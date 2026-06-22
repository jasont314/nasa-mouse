suppressPackageStartupMessages({
  library(DESeq2)
})

parse_args <- function(values) {
  parsed <- list()
  index <- 1
  while (index <= length(values)) {
    key <- values[[index]]
    if (!startsWith(key, "--") || index == length(values)) {
      stop(paste("Invalid argument sequence at", key))
    }
    parsed[[substring(key, 3)]] <- values[[index + 1]]
    index <- index + 2
  }
  parsed
}

bh_adjust <- function(p_values) {
  p.adjust(p_values, method = "BH")
}

read_required <- function(args, name) {
  if (!name %in% names(args)) {
    stop(paste("Missing required argument:", name))
  }
  args[[name]]
}

write_result_table <- function(result, gene_symbols, accession, output_path, alpha, lfc_cutoff) {
  table <- as.data.frame(result)
  table$gene_id <- rownames(table)
  table$gene_symbol <- gene_symbols[table$gene_id]
  table$accession <- accession
  table$significant_padj05_abs_lfc1 <- (
    !is.na(table$padj) &
      table$padj < alpha &
      abs(table$log2FoldChange) >= lfc_cutoff
  )
  table$direction <- ifelse(
    table$log2FoldChange > 0,
    "up",
    ifelse(table$log2FoldChange < 0, "down", "flat")
  )
  table <- table[, c(
    "gene_id", "gene_symbol", "accession", "baseMean", "log2FoldChange",
    "lfcSE", "stat", "pvalue", "padj", "significant_padj05_abs_lfc1",
    "direction"
  )]
  write.table(
    table,
    output_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  table
}

run_accession_deseq2 <- function(counts, metadata, gene_symbols, accession, output_dir, alpha, lfc_cutoff, min_count, min_samples) {
  selected_metadata <- metadata[metadata$accession == accession, , drop = FALSE]
  selected_metadata <- selected_metadata[selected_metadata$condition %in% c("FLT", "GC"), , drop = FALSE]
  condition_counts <- table(selected_metadata$condition)
  if (!all(c("FLT", "GC") %in% names(condition_counts)) ||
      condition_counts[["FLT"]] < 2 ||
      condition_counts[["GC"]] < 2) {
    message(sprintf("Skipping %s: fewer than 2 FLT or GC samples", accession))
    return(NULL)
  }

  selected_counts <- counts[, selected_metadata$sample, drop = FALSE]
  keep <- rowSums(selected_counts >= min_count) >= min_samples
  filtered_counts <- selected_counts[keep, , drop = FALSE]
  if (nrow(filtered_counts) == 0) {
    message(sprintf("Skipping %s: no genes passed count filter", accession))
    return(NULL)
  }

  rownames(selected_metadata) <- selected_metadata$sample
  selected_metadata$condition <- relevel(factor(selected_metadata$condition), ref = "GC")
  dds <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = selected_metadata,
    design = ~ condition
  )
  dds <- DESeq(dds, quiet = TRUE)
  result <- results(
    dds,
    contrast = c("condition", "FLT", "GC"),
    alpha = alpha
  )

  output_path <- file.path(output_dir, paste0(accession, "_flight_vs_ground.tsv"))
  table <- write_result_table(result, gene_symbols, accession, output_path, alpha, lfc_cutoff)
  table$n_flt <- as.integer(condition_counts[["FLT"]])
  table$n_gc <- as.integer(condition_counts[["GC"]])
  table$genes_retained_in_study <- nrow(filtered_counts)
  table
}

fixed_effect_meta <- function(per_study, alpha, lfc_cutoff, min_studies) {
  split_by_gene <- split(per_study, per_study$gene_id)
  rows <- lapply(split_by_gene, function(table) {
    valid <- (
      !is.na(table$log2FoldChange) &
        !is.na(table$lfcSE) &
        is.finite(table$log2FoldChange) &
        is.finite(table$lfcSE) &
        table$lfcSE > 0
    )
    tested <- table[valid, , drop = FALSE]
    if (nrow(tested) == 0) {
      return(NULL)
    }
    weights <- 1 / (tested$lfcSE^2)
    weight_sum <- sum(weights)
    meta_lfc <- sum(weights * tested$log2FoldChange) / weight_sum
    meta_se <- sqrt(1 / weight_sum)
    meta_z <- meta_lfc / meta_se
    meta_p <- 2 * pnorm(abs(meta_z), lower.tail = FALSE)
    q_stat <- sum(weights * ((tested$log2FoldChange - meta_lfc)^2))
    df <- nrow(tested) - 1
    i2 <- if (df > 0 && q_stat > 0) max(0, (q_stat - df) / q_stat) * 100 else NA_real_
    heterogeneity_p <- if (df > 0) pchisq(q_stat, df = df, lower.tail = FALSE) else NA_real_
    data.frame(
      gene_id = tested$gene_id[[1]],
      gene_symbol = tested$gene_symbol[[1]],
      studies_tested = nrow(tested),
      meta_log2_fold_change = meta_lfc,
      meta_se = meta_se,
      meta_z = meta_z,
      meta_p_value = meta_p,
      direction_up_studies = sum(tested$log2FoldChange > 0),
      direction_down_studies = sum(tested$log2FoldChange < 0),
      q_heterogeneity = q_stat,
      i2_percent = i2,
      heterogeneity_p_value = heterogeneity_p,
      stringsAsFactors = FALSE
    )
  })
  meta <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  meta$meta_fdr_bh <- bh_adjust(meta$meta_p_value)
  meta$eligible_meta <- meta$studies_tested >= min_studies
  meta$significant_fdr05_abs_log2fc1 <- (
    meta$eligible_meta &
    !is.na(meta$meta_fdr_bh) &
      meta$meta_fdr_bh < alpha &
      abs(meta$meta_log2_fold_change) >= lfc_cutoff
  )
  meta$consistent_direction <- ifelse(
    meta$direction_up_studies == meta$studies_tested,
    "all_up",
    ifelse(
      meta$direction_down_studies == meta$studies_tested,
      "all_down",
      "mixed"
    )
  )
  meta <- meta[order(meta$meta_fdr_bh, meta$meta_p_value), , drop = FALSE]
  rownames(meta) <- NULL
  meta
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
counts_path <- read_required(args, "counts")
metadata_path <- read_required(args, "metadata")
symbols_path <- read_required(args, "gene-symbols")
output_dir <- read_required(args, "output-dir")
alpha <- if ("alpha" %in% names(args)) as.numeric(args$alpha) else 0.05
lfc_cutoff <- if ("lfc-cutoff" %in% names(args)) as.numeric(args$`lfc-cutoff`) else 1
min_count <- if ("min-count" %in% names(args)) as.numeric(args$`min-count`) else 10
min_samples <- if ("min-samples" %in% names(args)) as.numeric(args$`min-samples`) else 3
min_studies <- if ("min-studies" %in% names(args)) as.numeric(args$`min-studies`) else 2

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
per_study_dir <- file.path(output_dir, "per_study")
dir.create(per_study_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading aggregate raw counts")
counts <- read.delim(
  counts_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

metadata <- read.delim(metadata_path, stringsAsFactors = FALSE)
symbols <- read.delim(symbols_path, stringsAsFactors = FALSE)
gene_symbols <- setNames(symbols$gene_symbol, symbols$gene_id)

missing_samples <- setdiff(metadata$sample, colnames(counts))
if (length(missing_samples) > 0) {
  stop(paste("Counts are missing samples:", paste(missing_samples, collapse = ", ")))
}

accessions <- sort(unique(metadata$accession))
message(sprintf("Running per-study DESeq2 for %d accessions", length(accessions)))
tables <- list()
for (accession in accessions) {
  message(sprintf("DESeq2 %s", accession))
  table <- run_accession_deseq2(
    counts,
    metadata,
    gene_symbols,
    accession,
    per_study_dir,
    alpha,
    lfc_cutoff,
    min_count,
    min_samples
  )
  if (!is.null(table)) {
    tables[[accession]] <- table
  }
}

if (length(tables) == 0) {
  stop("No accession-level DESeq2 result tables were produced")
}

per_study <- do.call(rbind, tables)
rownames(per_study) <- NULL
write.table(
  per_study,
  file.path(output_dir, "per_study_deseq2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Combining per-study DESeq2 effects")
meta <- fixed_effect_meta(per_study, alpha, lfc_cutoff, min_studies)
write.table(
  meta,
  file.path(output_dir, "deseq2_meta.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  head(meta[meta$eligible_meta, , drop = FALSE], 200),
  file.path(output_dir, "top_deseq2_meta_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary <- data.frame(
  accessions_tested = paste(names(tables), collapse = ","),
  genes_tested = nrow(meta),
  genes_eligible_min_studies = sum(meta$eligible_meta, na.rm = TRUE),
  min_studies = min_studies,
  significant_fdr05_abs_log2fc1 = sum(meta$significant_fdr05_abs_log2fc1, na.rm = TRUE),
  significant_up = sum(
    meta$significant_fdr05_abs_log2fc1 &
      meta$meta_log2_fold_change > 0,
    na.rm = TRUE
  ),
  significant_down = sum(
    meta$significant_fdr05_abs_log2fc1 &
      meta$meta_log2_fold_change < 0,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)
write.table(
  summary,
  file.path(output_dir, "deseq2_meta_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top_sig <- meta[meta$significant_fdr05_abs_log2fc1, , drop = FALSE]
top_sig <- top_sig[order(-abs(top_sig$meta_log2_fold_change)), , drop = FALSE]
top_sig <- head(top_sig[, c(
  "gene_symbol", "gene_id", "studies_tested", "meta_log2_fold_change",
  "meta_fdr_bh", "direction_up_studies", "direction_down_studies",
  "consistent_direction", "i2_percent"
)], 15)
top_sig_lines <- capture.output(write.table(
  top_sig,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
))
report <- c(
  "# Aggregate Liver DESeq2 Meta-Analysis",
  "",
  "Per-study DESeq2 was run separately for each OSD accession, then",
  "log2 fold-change estimates were combined with a fixed-effect inverse-variance",
  "meta-analysis. The aggregate significance flag requires the gene to be tested",
  sprintf("in at least %d accessions.", min_studies),
  "",
  sprintf("- Accessions tested: %s", paste(names(tables), collapse = ", ")),
  sprintf("- Genes tested in at least one study: %s", format(nrow(meta), big.mark = ",")),
  sprintf(
    "- Genes eligible for aggregate meta-analysis: %s",
    format(sum(meta$eligible_meta, na.rm = TRUE), big.mark = ",")
  ),
  sprintf(
    "- Significant genes at FDR < %.2f and abs(log2FC) >= %.1f: %d",
    alpha,
    lfc_cutoff,
    sum(meta$significant_fdr05_abs_log2fc1, na.rm = TRUE)
  ),
  sprintf(
    "- Up: %d",
    sum(
      meta$significant_fdr05_abs_log2fc1 &
        meta$meta_log2_fold_change > 0,
      na.rm = TRUE
    )
  ),
  sprintf(
    "- Down: %d",
    sum(
      meta$significant_fdr05_abs_log2fc1 &
        meta$meta_log2_fold_change < 0,
      na.rm = TRUE
    )
  ),
  "",
  "## Largest Significant Effects",
  "",
  "```tsv",
  top_sig_lines,
  "```"
)
writeLines(report, file.path(output_dir, "DESEQ2_META_SUMMARY.md"))
message("Done")
