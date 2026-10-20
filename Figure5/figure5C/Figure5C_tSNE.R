#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required_packages <- c("Rtsne", "ggplot2", "vegan", "permute")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Rtsne)
  library(ggplot2)
  library(vegan)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) return(dirname(normalizePath(active_path)))
  }
  getwd()
}

parse_options <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) stop("Arguments must be supplied as --name value pairs.")
  out <- list()
  if (length(args) == 0) return(out)
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    if (!nzchar(key) || grepl("^--", args[i + 1])) stop("Invalid command-line arguments.")
    out[[key]] <- args[i + 1]
  }
  out
}

cli <- parse_options()
base_dir <- get_script_dir()
matrix_file <- if (!is.null(cli$matrix)) cli$matrix else
  file.path(base_dir, "virus_table_selected_samples_updated.csv")
metadata_file <- if (!is.null(cli$metadata)) cli$metadata else
  file.path(base_dir, "library_annotation.csv")
output_dir <- if (!is.null(cli$outdir)) cli$outdir else base_dir
requested_perplexity <- if (!is.null(cli$perplexity)) as.numeric(cli$perplexity) else 7
seed <- if (!is.null(cli$seed)) as.integer(cli$seed) else 123
max_iter <- if (!is.null(cli$max_iter)) as.integer(cli$max_iter) else 3000
n_permutations <- if (!is.null(cli$permutations)) as.integer(cli$permutations) else 999L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(matrix_file)) stop("Missing matrix input: ", matrix_file)
if (!file.exists(metadata_file)) stop("Missing metadata input: ", metadata_file)
if (!is.finite(requested_perplexity) || requested_perplexity <= 0) {
  stop("Perplexity must be a positive number.")
}
if (!is.finite(seed)) stop("Seed must be an integer.")
if (!is.finite(max_iter) || max_iter < 500) stop("max_iter must be at least 500.")
if (!is.finite(n_permutations) || n_permutations < 99) {
  stop("permutations must be an integer of at least 99.")
}

matrix_raw <- read.csv(matrix_file, check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
metadata_raw <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA"))

if (!"species_name" %in% names(matrix_raw)) stop("The matrix must contain species_name.")
metadata_required <- c("lib_id", "Location", "Domestic or Wild", "Sample Type", "Order")
missing_metadata <- setdiff(metadata_required, names(metadata_raw))
if (length(missing_metadata) > 0) {
  stop("Metadata is missing: ", paste(missing_metadata, collapse = ", "))
}

sample_columns <- setdiff(names(matrix_raw), "species_name")
if (length(sample_columns) < 3) stop("The matrix contains too few library columns.")
if (anyDuplicated(metadata_raw$lib_id)) stop("Metadata contains duplicated lib_id values.")
missing_in_metadata <- setdiff(sample_columns, metadata_raw$lib_id)
if (length(missing_in_metadata) > 0) {
  stop("Libraries missing from metadata: ", paste(missing_in_metadata, collapse = ", "))
}

species_labels <- trimws(as.character(matrix_raw$species_name))
if (anyNA(species_labels) || any(!nzchar(species_labels))) {
  stop("species_name contains missing or empty values.")
}

matrix_values <- as.matrix(matrix_raw[sample_columns])
storage.mode(matrix_values) <- "numeric"
matrix_values[is.na(matrix_values)] <- 0
matrix_values[matrix_values < 0] <- 0
aggregated_values <- rowsum(matrix_values, group = species_labels, reorder = FALSE)
matrix_aggregated <- data.frame(
  species_name = rownames(aggregated_values),
  aggregated_values,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(matrix_aggregated) <- NULL

metadata <- metadata_raw[match(sample_columns, metadata_raw$lib_id), metadata_required]
rownames(metadata) <- metadata$lib_id

feature_matrix <- do.call(cbind, lapply(matrix_aggregated[sample_columns], function(x) {
  as.numeric(as.character(x))
}))
colnames(feature_matrix) <- sample_columns
rownames(feature_matrix) <- matrix_aggregated$species_name

mixed_mask <- !is.na(metadata$`Sample Type`) & metadata$`Sample Type` == "Mixed type libs"
if (any(mixed_mask)) {
  feature_matrix <- feature_matrix[, !mixed_mask, drop = FALSE]
  metadata <- metadata[!mixed_mask, , drop = FALSE]
}

raw_sample_matrix <- t(feature_matrix)
raw_sample_matrix[is.na(raw_sample_matrix)] <- 0
raw_sample_matrix[raw_sample_matrix < 0] <- 0

zero_library_mask <- rowSums(raw_sample_matrix) == 0
zero_library_ids <- rownames(raw_sample_matrix)[zero_library_mask]
if (length(zero_library_ids) > 0) {
  raw_sample_matrix <- raw_sample_matrix[!zero_library_mask, , drop = FALSE]
  metadata <- metadata[!metadata$lib_id %in% zero_library_ids, , drop = FALSE]
}

zero_features <- colSums(raw_sample_matrix) == 0
raw_sample_matrix <- raw_sample_matrix[, !zero_features, drop = FALSE]

reference_sample_types <- c(
  "Blood", "Brain", "Cecal content", "Cloacal swab", "Fecal samples",
  "Intestine", "Lesion", "Liver", "Lung"
)
sample_type_values <- trimws(as.character(metadata$`Sample Type`))
keep_sample_type <- sample_type_values %in% reference_sample_types
if (any(!keep_sample_type)) {
  raw_sample_matrix <- raw_sample_matrix[keep_sample_type, , drop = FALSE]
  metadata <- metadata[keep_sample_type, , drop = FALSE]
  sample_type_values <- sample_type_values[keep_sample_type]
}

observed_sample_types <- unique(sample_type_values)
sample_type_levels <- reference_sample_types[reference_sample_types %in% observed_sample_types]
metadata$Sample_Type <- factor(sample_type_values, levels = sample_type_levels)
metadata$Block_Order <- factor(trimws(as.character(metadata$Order)))
metadata$Block_Location <- factor(trimws(as.character(metadata$Location)))
rownames(metadata) <- metadata$lib_id

if (nrow(raw_sample_matrix) < 8) stop("Too few libraries remain for t-SNE.")
if (ncol(raw_sample_matrix) < 2) stop("Too few non-zero features remain for t-SNE.")

sample_ids <- rownames(raw_sample_matrix)
metadata <- metadata[match(sample_ids, metadata$lib_id), , drop = FALSE]
log_sample_matrix <- log10(raw_sample_matrix + 1)

bray_curtis <- vegan::vegdist(log_sample_matrix, method = "bray")

make_permutation_control <- function(data, block_name = NULL) {
  if (is.null(block_name)) return(n_permutations)
  block_values <- droplevels(factor(data[[block_name]]))
  permute::how(nperm = n_permutations, blocks = block_values)
}

run_permanova <- function(distance, data, block_name = NULL) {
  block_label <- if (is.null(block_name)) "unblocked" else block_name
  result <- data.frame(
    test = "PERMANOVA", block = block_label, term = "Sample_Type",
    n = nrow(data), n_blocks = if (is.null(block_name)) NA_integer_ else length(unique(data[[block_name]])),
    permutations = n_permutations,
    df = NA_real_, R2 = NA_real_, F = NA_real_, p_value = NA_real_,
    status = "failed", stringsAsFactors = FALSE
  )
  if (nlevels(droplevels(data$Sample_Type)) < 2) {
    result$status <- "fewer than two sample types"
    return(result)
  }
  fit <- tryCatch({
    vegan::adonis2(
      distance ~ Sample_Type,
      data = data,
      permutations = make_permutation_control(data, block_name),
      by = "terms"
    )
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    result$status <- conditionMessage(fit)
    return(result)
  }
  term_row <- fit["Sample_Type", , drop = FALSE]
  result$df <- unname(term_row$Df)
  result$R2 <- unname(term_row$R2)
  result$F <- unname(term_row$F)
  result$p_value <- unname(term_row$`Pr(>F)`)
  result$status <- "ok"
  result
}

permanova_results <- rbind(
  run_permanova(bray_curtis, metadata),
  run_permanova(bray_curtis, metadata, "Block_Order"),
  run_permanova(bray_curtis, metadata, "Block_Location")
)

run_permdisp <- function(distance, data, block_name = NULL) {
  block_label <- if (is.null(block_name)) "unblocked" else block_name
  result <- data.frame(
    test = "PERMDISP", block = block_label, grouping = "Sample_Type",
    n = nrow(data), n_blocks = if (is.null(block_name)) NA_integer_ else length(unique(data[[block_name]])),
    permutations = n_permutations,
    F = NA_real_, p_value = NA_real_, status = "failed", stringsAsFactors = FALSE
  )
  if (nlevels(droplevels(data$Sample_Type)) < 2) {
    result$status <- "fewer than two sample types"
    return(result)
  }
  dispersion <- tryCatch(
    vegan::betadisper(distance, data$Sample_Type),
    error = function(e) e
  )
  if (inherits(dispersion, "error")) {
    result$status <- conditionMessage(dispersion)
    return(result)
  }
  dispersion_test <- tryCatch(
    vegan::permutest(
      dispersion,
      permutations = make_permutation_control(data, block_name)
    ),
    error = function(e) e
  )
  if (inherits(dispersion_test, "error")) {
    result$status <- conditionMessage(dispersion_test)
    return(result)
  }
  result$F <- unname(dispersion_test$tab[1, "F"])
  result$p_value <- unname(dispersion_test$tab[1, "Pr(>F)"])
  result$status <- "ok"
  result
}

permdisp_results <- rbind(
  run_permdisp(bray_curtis, metadata),
  run_permdisp(bray_curtis, metadata, "Block_Order"),
  run_permdisp(bray_curtis, metadata, "Block_Location")
)

max_allowed_perplexity <- floor((nrow(log_sample_matrix) - 1) / 3)
if (max_allowed_perplexity < 2) stop("Too few libraries for the requested t-SNE perplexity.")
perplexity <- min(requested_perplexity, max_allowed_perplexity)

set.seed(seed)
tsne_fit <- Rtsne::Rtsne(
  X = log_sample_matrix,
  dims = 2,
  perplexity = perplexity,
  theta = 0.5,
  max_iter = max_iter,
  pca = TRUE,
  check_duplicates = FALSE,
  verbose = FALSE
)

# t-SNE does not have the NMDS stress statistic.  Rtsne reports the final
# Kullback–Leibler objective instead; this is the appropriate optimization
# diagnostic to report for a t-SNE embedding.
final_kl <- if (!is.null(tsne_fit$itercosts) && length(tsne_fit$itercosts) > 0) {
  as.numeric(tail(tsne_fit$itercosts, 1))
} else {
  sum(as.numeric(tsne_fit$costs), na.rm = TRUE)
}
if (!is.finite(final_kl)) final_kl <- NA_real_
final_kl_per_library <- if (is.finite(final_kl)) final_kl / nrow(log_sample_matrix) else NA_real_

tsne_coordinates <- data.frame(
  lib_id = sample_ids,
  tSNE1 = tsne_fit$Y[, 1],
  tSNE2 = tsne_fit$Y[, 2],
  metadata[match(sample_ids, metadata$lib_id),
           c("Sample_Type", "Location", "Domestic or Wild", "Order"), drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_palette <- c(
  Blood = "#DD5F60",
  Brain = "#238B45",
  "Cecal content" = "#6A51A3",
  "Cloacal swab" = "#EC7014",
  "Fecal samples" = "#3690C0",
  Intestine = "#DD3497",
  Lesion = "#7BCCC4",
  Liver = "#C79494",
  Lung = "#D95F02"
)
sample_palette <- sample_palette[sample_type_levels]

format_p <- function(value) {
  if (is.na(value)) return("NA")
  if (value < 0.001) return("<0.001")
  formatC(value, format = "f", digits = 3)
}

blocked_order <- permanova_results[permanova_results$block == "Block_Order", , drop = FALSE]
blocked_order_disp <- permdisp_results[permdisp_results$block == "Block_Order", , drop = FALSE]
stat_label <- paste0(
  "PERMANOVA: R2 = ",
  ifelse(nrow(blocked_order) == 1, formatC(blocked_order$R2, format = "f", digits = 3), "NA"),
  ", P = ", ifelse(nrow(blocked_order) == 1, format_p(blocked_order$p_value), "NA"),
  "\nPERMDISP: P = ",
  ifelse(nrow(blocked_order_disp) == 1, format_p(blocked_order_disp$p_value), "NA"),
  "\nt-SNE KL = ", ifelse(is.finite(final_kl), formatC(final_kl, format = "f", digits = 3), "NA")
)

tsne_plot <- ggplot(tsne_coordinates,
                    aes(x = tSNE1, y = tSNE2, colour = Sample_Type)) +
  geom_point(shape = 16, size = 4.6, alpha = 0.52) +
  annotate("text", x = -Inf, y = Inf,
           label = stat_label,
           hjust = -0.05, vjust = 1.15, size = 5.2, colour = "black") +
  scale_colour_manual(values = sample_palette, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 5.8, alpha = 1))) +
  labs(x = "t-SNE1", y = "t-SNE2", colour = "Sample type") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.key.width = grid::unit(1.25, "cm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(8, 12, 8, 72)
  )

write.csv(tsne_coordinates, file.path(output_dir, "tsne_coordinates.csv"), row.names = FALSE)
write.csv(metadata, file.path(output_dir, "tsne_metadata_used.csv"), row.names = FALSE)
write.csv(permanova_results, file.path(output_dir, "tsne_permanova_results.csv"), row.names = FALSE)
write.csv(permdisp_results, file.path(output_dir, "tsne_permdisp_results.csv"), row.names = FALSE)
write.csv(data.frame(
  metric = c("final_tSNE_KL_divergence", "final_tSNE_KL_divergence_per_library"),
  value = c(final_kl, final_kl_per_library),
  interpretation = c(
    "Rtsne final optimization objective; not NMDS stress",
    "Final t-SNE KL divergence divided by the number of libraries"
  ),
  stringsAsFactors = FALSE
), file.path(output_dir, "tsne_optimization_metrics.csv"), row.names = FALSE)
write.csv(data.frame(lib_id = zero_library_ids,
                     reason = "zero abundance across the input matrix"),
          file.path(output_dir, "zero_abundance_libraries_excluded.csv"), row.names = FALSE)

summary_lines <- c(
  paste0("Input matrix: ", basename(matrix_file)),
  paste0("Metadata: ", basename(metadata_file)),
  paste0("Libraries analysed: ", nrow(tsne_coordinates)),
  paste0("Features analysed: ", ncol(raw_sample_matrix)),
  paste0("Merged duplicate species rows: ", sum(table(species_labels) > 1)),
  paste0("Mixed-tissue libraries excluded: ", sum(mixed_mask)),
  paste0("Zero-abundance libraries excluded: ", length(zero_library_ids)),
  paste0("Sample types analysed: ", paste(sample_type_levels, collapse = "; ")),
  "Transformation: log10(RPM+1)",
  paste0("Perplexity requested: ", requested_perplexity),
  paste0("Perplexity used: ", perplexity),
  paste0("Seed: ", seed),
  paste0("Max iterations: ", max_iter),
  paste0("Permutation count: ", n_permutations),
  "Blocked permutations: permute::how(blocks = Order) and permute::how(blocks = Location)",
  paste0("Order-stratified PERMANOVA R2: ",
         ifelse(nrow(blocked_order) == 1, formatC(blocked_order$R2, format = "f", digits = 6), "NA")),
  paste0("Order-stratified PERMANOVA P: ",
         ifelse(nrow(blocked_order) == 1, format_p(blocked_order$p_value), "NA")),
  paste0("Order-blocked PERMDISP P: ",
         ifelse(nrow(blocked_order_disp) == 1, format_p(blocked_order_disp$p_value), "NA")),
  paste0("t-SNE final KL divergence: ",
         ifelse(is.finite(final_kl), formatC(final_kl, format = "f", digits = 6), "NA")),
  "Independent orthogonal data and a pre-specified biological negative-control variable were not supplied; no such result is claimed.",
  "t-SNE is exploratory and does not replace PERMANOVA or PERMDISP."
)
writeLines(summary_lines, file.path(output_dir, "tsne_analysis_summary.txt"))

ggsave(file.path(output_dir, "Figure5C_tSNE.pdf"), tsne_plot,
       width = 10, height = 7, units = "in", device = "pdf")
ggsave(file.path(output_dir, "Figure5C_tSNE.png"), tsne_plot,
       width = 10, height = 7, units = "in", dpi = 300)

capture.output(sessionInfo(), file = file.path(output_dir, "tsne_sessionInfo.txt"))
message("t-SNE completed: ", output_dir)
