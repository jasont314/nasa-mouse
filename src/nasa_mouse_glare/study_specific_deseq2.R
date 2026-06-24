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

read_required <- function(args, name) {
  if (!name %in% names(args)) {
    stop(paste("Missing required argument:", name))
  }
  args[[name]]
}

run_deseq_with_dispersion_fallback <- function(dds, accession) {
  tryCatch(
    {
      dds <- DESeq(dds, quiet = TRUE)
      metadata(dds)$dispersion_fit <- "default"
      dds
    },
    error = function(error) {
      message_text <- conditionMessage(error)
      if (!grepl("all gene-wise dispersion estimates", message_text, fixed = TRUE)) {
        stop(error)
      }
      message(sprintf(
        "DESeq2 %s: default dispersion fit failed; using gene-wise dispersion estimates",
        accession
      ))
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersionsGeneEst(dds, quiet = TRUE)
      dispersions(dds) <- mcols(dds)$dispGeneEst
      dds <- nbinomWaldTest(dds, quiet = TRUE)
      metadata(dds)$dispersion_fit <- "gene_wise_fallback"
      dds
    }
  )
}

write_result_table <- function(result, gene_symbols, accession, output_path, alpha) {
  table <- as.data.frame(result)
  table$gene_id <- rownames(table)
  table$gene_symbol <- gene_symbols[table$gene_id]
  table$accession <- accession
  table$significant_padj05 <- !is.na(table$padj) & table$padj < alpha
  table$direction <- ifelse(
    table$log2FoldChange > 0,
    "up_in_flight",
    ifelse(table$log2FoldChange < 0, "down_in_flight", "flat")
  )
  table <- table[, c(
    "gene_id", "gene_symbol", "accession", "baseMean", "log2FoldChange",
    "lfcSE", "stat", "pvalue", "padj", "significant_padj05", "direction"
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

choose_design <- function(selected_metadata) {
  selected_metadata$stratum <- factor(selected_metadata$stratum)
  strata <- levels(selected_metadata$stratum)
  if (length(strata) < 2) {
    return(list(formula = ~ condition, name = "~ condition"))
  }

  condition_by_stratum <- table(selected_metadata$stratum, selected_metadata$condition)
  balanced <- all(condition_by_stratum[, "flight"] > 0) &&
    all(condition_by_stratum[, "ground"] > 0)
  if (!balanced) {
    return(list(formula = ~ condition, name = "~ condition"))
  }

  design <- model.matrix(~ stratum + condition, data = selected_metadata)
  if (qr(design)$rank < ncol(design)) {
    return(list(formula = ~ condition, name = "~ condition"))
  }
  list(formula = ~ stratum + condition, name = "~ stratum + condition")
}

run_accession <- function(counts, metadata, gene_symbols, accession, output_dir, alpha, min_count, min_samples) {
  selected_metadata <- metadata[metadata$accession == accession, , drop = FALSE]
  selected_metadata <- selected_metadata[selected_metadata$condition %in% c("flight", "ground"), , drop = FALSE]
  condition_counts <- table(selected_metadata$condition)
  if (!all(c("flight", "ground") %in% names(condition_counts)) ||
      condition_counts[["flight"]] < 2 ||
      condition_counts[["ground"]] < 2) {
    message(sprintf("Skipping %s: fewer than 2 flight or ground samples", accession))
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
  selected_metadata$condition <- relevel(factor(selected_metadata$condition), ref = "ground")
  selected_metadata$stratum <- factor(selected_metadata$stratum)
  chosen_design <- choose_design(selected_metadata)

  dds <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = selected_metadata,
    design = chosen_design$formula
  )
  dds <- run_deseq_with_dispersion_fallback(dds, accession)
  result <- results(
    dds,
    contrast = c("condition", "flight", "ground"),
    alpha = alpha
  )

  output_path <- file.path(output_dir, paste0(accession, "_flight_vs_ground.tsv"))
  table <- write_result_table(result, gene_symbols, accession, output_path, alpha)
  table$n_flight <- as.integer(condition_counts[["flight"]])
  table$n_ground <- as.integer(condition_counts[["ground"]])
  table$genes_retained_in_study <- nrow(filtered_counts)
  table$design <- chosen_design$name
  table$dispersion_fit <- metadata(dds)$dispersion_fit
  table
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
counts_path <- read_required(args, "counts")
metadata_path <- read_required(args, "metadata")
symbols_path <- read_required(args, "gene-symbols")
output_dir <- read_required(args, "output-dir")
alpha <- if ("alpha" %in% names(args)) as.numeric(args$alpha) else 0.05
min_count <- if ("min-count" %in% names(args)) as.numeric(args$`min-count`) else 10
min_samples <- if ("min-samples" %in% names(args)) as.numeric(args$`min-samples`) else 3

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
per_study_dir <- file.path(output_dir, "per_study")
dir.create(per_study_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading raw counts")
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
  table <- run_accession(
    counts,
    metadata,
    gene_symbols,
    accession,
    per_study_dir,
    alpha,
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

summary_rows <- lapply(names(tables), function(accession) {
  table <- tables[[accession]]
  data.frame(
    accession = accession,
    n_flight = table$n_flight[[1]],
    n_ground = table$n_ground[[1]],
    genes_tested = nrow(table),
    significant_padj05 = sum(table$significant_padj05, na.rm = TRUE),
    significant_up = sum(table$significant_padj05 & table$log2FoldChange > 0, na.rm = TRUE),
    significant_down = sum(table$significant_padj05 & table$log2FoldChange < 0, na.rm = TRUE),
    design = table$design[[1]],
    dispersion_fit = table$dispersion_fit[[1]],
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
write.table(
  summary,
  file.path(output_dir, "study_deseq2_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Done")
