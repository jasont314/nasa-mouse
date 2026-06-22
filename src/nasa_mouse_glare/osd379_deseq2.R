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

write_results <- function(result, output_path, alpha, lfc_cutoff) {
  table <- as.data.frame(result)
  table$gene_id <- rownames(table)
  table$significant_padj05_abs_lfc1 <- (
    !is.na(table$padj) &
      table$padj < alpha &
      abs(table$log2FoldChange) >= lfc_cutoff
  )
  table$direction <- ifelse(
    table$significant_padj05_abs_lfc1,
    ifelse(table$log2FoldChange > 0, "up", "down"),
    "not_deg"
  )
  table <- table[, c(
    "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
    "pvalue", "padj", "significant_padj05_abs_lfc1", "direction"
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

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("raw-counts", "matched-slots", "exclude-samples", "output-dir")
missing <- required[!required %in% names(args)]
if (length(missing) > 0) {
  stop(paste("Missing required arguments:", paste(missing, collapse = ", ")))
}

alpha <- if ("alpha" %in% names(args)) as.numeric(args$alpha) else 0.05
lfc_cutoff <- if ("lfc-cutoff" %in% names(args)) {
  as.numeric(args$`lfc-cutoff`)
} else {
  1
}
output_dir <- args$`output-dir`
filter_mode <- if ("filter-mode" %in% names(args)) args$`filter-mode` else "independent"
if (!filter_mode %in% c("independent", "matched")) {
  stop("filter-mode must be independent or matched")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading raw RSEM counts")
raw <- read.csv(
  args$`raw-counts`,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
raw <- as.matrix(raw)
storage.mode(raw) <- "integer"

slots <- read.delim(args$`matched-slots`, stringsAsFactors = FALSE)
exclusions <- read.delim(args$`exclude-samples`, stringsAsFactors = FALSE)
if ("sample" %in% colnames(exclusions)) {
  excluded <- exclusions$sample
} else if ("directly_flagged_profile" %in% colnames(exclusions)) {
  excluded <- exclusions$directly_flagged_profile
} else {
  stop("Exclusion table requires sample or directly_flagged_profile")
}
excluded <- unique(excluded[!is.na(excluded) & excluded != ""])
drop_slot <- slots$flt_profile %in% excluded | slots$gc_profile %in% excluded
removed_slots <- slots[drop_slot, , drop = FALSE]
if (filter_mode == "matched") {
  analysis_slots <- slots[!drop_slot, , drop = FALSE]
} else {
  analysis_slots <- slots
}
if (nrow(analysis_slots) == 0) {
  stop("No matched slots remain after filtering")
}
write.table(
  removed_slots,
  file.path(output_dir, "excluded_matched_feature_slots.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

make_metadata <- function(slots, excluded, filter_mode) {
  flight <- data.frame(
    sample = slots$flt_profile,
    condition = "flight",
    stratum = paste(slots$cohort, slots$age, sep = "_"),
    feature = slots$feature,
    stringsAsFactors = FALSE
  )
  ground <- data.frame(
    sample = slots$gc_profile,
    condition = "ground",
    stratum = paste(slots$cohort, slots$age, sep = "_"),
    feature = slots$feature,
    stringsAsFactors = FALSE
  )
  if (filter_mode == "independent") {
    flight <- flight[!flight$sample %in% excluded, , drop = FALSE]
    ground <- ground[!ground$sample %in% excluded, , drop = FALSE]
  }
  metadata <- rbind(flight, ground)
  metadata$condition <- relevel(factor(metadata$condition), ref = "ground")
  metadata$stratum <- factor(metadata$stratum)
  rownames(metadata) <- metadata$sample
  metadata
}

metadata <- make_metadata(analysis_slots, excluded, filter_mode)
use_composition_covariate <- "composition-scores" %in% names(args)
if (use_composition_covariate) {
  composition <- read.delim(
    args$`composition-scores`,
    stringsAsFactors = FALSE
  )
  required_composition <- c("sample", "muscle_mean_log2")
  missing_composition <- setdiff(required_composition, colnames(composition))
  if (length(missing_composition) > 0) {
    stop(paste(
      "Composition table is missing:",
      paste(missing_composition, collapse = ", ")
    ))
  }
  score_map <- setNames(composition$muscle_mean_log2, composition$sample)
  metadata$muscle_score <- score_map[metadata$sample]
  if (any(is.na(metadata$muscle_score))) {
    stop("Composition score is missing for one or more retained samples")
  }
  metadata$muscle_score_z <- ave(
    metadata$muscle_score,
    metadata$stratum,
    FUN = function(values) as.numeric(scale(values))
  )
}
missing_samples <- setdiff(metadata$sample, colnames(raw))
if (length(missing_samples) > 0) {
  stop(paste("Raw-count table is missing samples:", paste(missing_samples, collapse = ", ")))
}
counts_all <- raw[, metadata$sample, drop = FALSE]
keep <- rowSums(counts_all >= 10) >= 3
counts <- counts_all[keep, , drop = FALSE]
message(
  sprintf(
    "Retained %d/%d genes and %d samples after filtering",
    nrow(counts), nrow(counts_all), ncol(counts)
  )
)

global_design <- if (use_composition_covariate) {
  ~ stratum + muscle_score_z + condition
} else {
  ~ stratum + condition
}
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = global_design
)
dds <- DESeq(dds)
global <- results(
  dds,
  contrast = c("condition", "flight", "ground"),
  alpha = alpha
)
global_table <- write_results(
  global,
  file.path(output_dir, "global_flight_vs_ground.tsv"),
  alpha,
  lfc_cutoff
)

size_factors <- sizeFactors(dds)
normalized_all <- sweep(counts_all, 2, size_factors, "/")
write.csv(
  normalized_all,
  file.path(output_dir, "filtered_deseq2_normalized_counts.csv"),
  quote = FALSE
)
if ("glare-gene-reference" %in% names(args)) {
  glare_reference <- read.csv(
    args$`glare-gene-reference`,
    row.names = 1,
    check.names = FALSE,
    nrows = nrow(raw),
    stringsAsFactors = FALSE
  )
  glare_genes <- intersect(rownames(normalized_all), rownames(glare_reference))
  write.csv(
    normalized_all[glare_genes, , drop = FALSE],
    file.path(output_dir, "filtered_deseq2_normalized_counts_glare.csv"),
    quote = FALSE
  )
}
write.table(
  data.frame(sample = names(size_factors), size_factor = size_factors),
  file.path(output_dir, "size_factors.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  as.data.frame(metadata),
  file.path(output_dir, "sample_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

vsd <- vst(dds, blind = FALSE)
pca <- prcomp(t(assay(vsd)))
variance <- pca$sdev^2 / sum(pca$sdev^2)
pca_table <- data.frame(
  sample = rownames(pca$x),
  condition = metadata[rownames(pca$x), "condition"],
  stratum = metadata[rownames(pca$x), "stratum"],
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC1_variance = variance[[1]],
  PC2_variance = variance[[2]]
)
write.table(
  pca_table,
  file.path(output_dir, "vst_pca.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

stratum_tables <- list()
for (stratum_name in levels(metadata$stratum)) {
  selected <- metadata$stratum == stratum_name
  stratum_metadata <- droplevels(metadata[selected, , drop = FALSE])
  stratum_counts_all <- counts_all[, selected, drop = FALSE]
  stratum_keep <- rowSums(stratum_counts_all >= 10) >= 3
  stratum_counts <- stratum_counts_all[stratum_keep, , drop = FALSE]
  stratum_design <- if (use_composition_covariate) {
    ~ muscle_score_z + condition
  } else {
    ~ condition
  }
  stratum_dds <- DESeqDataSetFromMatrix(
    countData = stratum_counts,
    colData = stratum_metadata,
    design = stratum_design
  )
  stratum_dds <- DESeq(stratum_dds)
  stratum_result <- results(
    stratum_dds,
    contrast = c("condition", "flight", "ground"),
    alpha = alpha
  )
  output_path <- file.path(
    output_dir,
    paste0("stratum_", gsub("[^A-Za-z0-9]+", "_", stratum_name), ".tsv")
  )
  table <- write_results(stratum_result, output_path, alpha, lfc_cutoff)
  table$stratum <- stratum_name
  stratum_tables[[stratum_name]] <- table
}

stratum_long <- do.call(rbind, stratum_tables)
write.table(
  stratum_long,
  file.path(output_dir, "stratum_results_long.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary <- data.frame(
  analysis = c("global_adjusted", names(stratum_tables)),
  samples = c(
    nrow(metadata),
    vapply(
      names(stratum_tables),
      function(name) sum(metadata$stratum == name),
      numeric(1)
    )
  ),
  genes_tested = c(
    nrow(global_table),
    vapply(stratum_tables, nrow, numeric(1))
  ),
  significant = c(
    sum(global_table$significant_padj05_abs_lfc1),
    vapply(
      stratum_tables,
      function(table) sum(table$significant_padj05_abs_lfc1),
      numeric(1)
    )
  ),
  up = c(
    sum(global_table$direction == "up"),
    vapply(stratum_tables, function(table) sum(table$direction == "up"), numeric(1))
  ),
  down = c(
    sum(global_table$direction == "down"),
    vapply(stratum_tables, function(table) sum(table$direction == "down"), numeric(1))
  )
)
summary$model <- c(
  paste(deparse(global_design), collapse = " "),
  rep(paste(deparse(stratum_design), collapse = " "), length(stratum_tables))
)
write.table(
  summary,
  file.path(output_dir, "deseq2_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message("DESeq2 analysis complete: ", output_dir)
