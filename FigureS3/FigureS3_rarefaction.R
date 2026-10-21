# Figure S3 rarefaction and diversity analysis.
# Inputs: FigureS3_libs_annot.xlsx and FigureS3_virus_table_filtered.csv.

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install it with install.packages('readxl').")
}

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (!is.null(script_file) && nzchar(script_file)) {
  dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
input_dir <- script_dir
annotation_file <- file.path(input_dir, "FigureS3_libs_annot.xlsx")
abundance_file <- file.path(input_dir, "FigureS3_virus_table_filtered.csv")

SEED <- 20260902L
N_ITER <- 10000L
TYPE_DEPTH <- 5L

TYPE_CLASS <- c(
  "Fecal samples" = "Mucosal/luminal",
  "Cloacal swab" = "Mucosal/luminal",
  "Cecal content" = "Mucosal/luminal",
  "Oral swab" = "Mucosal/luminal",
  "Liver" = "Solid tissues",
  "Intestine" = "Solid tissues",
  "Brain" = "Solid tissues",
  "Lung" = "Solid tissues",
  "Spleen" = "Solid tissues",
  "Kidney" = "Solid tissues",
  "Heart" = "Solid tissues",
  "Blood" = "Other",
  "Lesion" = "Other",
  "Feathers" = "Other"
)
TYPE_ORDER <- c(
  "Fecal samples", "Cloacal swab", "Cecal content", "Blood", "Lesion",
  "Liver", "Intestine", "Brain", "Lung"
)
CLASS_COLORS <- c(
  "Mucosal/luminal" = "#168A83",
  "Solid tissues" = "#D66B42",
  "Other" = "#687386"
)

annotation <- as.data.frame(
  readxl::read_excel(annotation_file, sheet = 1, .name_repair = "unique"),
  stringsAsFactors = FALSE
)
names(annotation) <- trimws(names(annotation))
sample_type_col <- "Sample Type"
if (!all(c("lib_id", sample_type_col) %in% names(annotation))) {
  stop("Annotation must contain 'lib_id' and 'Sample Type'.")
}
annotation$lib_id <- as.character(annotation$lib_id)
annotation[[sample_type_col]] <- trimws(as.character(annotation[[sample_type_col]]))

abundance_raw <- read.csv(
  abundance_file, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN")
)
metadata_columns <- c("lib_id", "kindom", "phylum", "family", "species_name")
if (!all(metadata_columns %in% names(abundance_raw))) {
  stop("Abundance matrix is missing one or more metadata columns.")
}
library_columns <- setdiff(names(abundance_raw), metadata_columns)
if (!setequal(library_columns, annotation$lib_id)) {
  stop("Library identifiers do not match between the annotation workbook and matrix.")
}

abundance <- as.data.frame(
  lapply(abundance_raw[library_columns], function(x) {
    z <- suppressWarnings(as.numeric(x))
    if (anyNA(z) || any(z < 0)) stop("Abundance matrix contains NA or negative values.")
    z
  }),
  check.names = FALSE, stringsAsFactors = FALSE
)
abundance <- as.matrix(abundance)
rownames(abundance) <- make.unique(as.character(abundance_raw$species_name))

false_positive_cells <- data.frame(
  species_name = c("Megrivirus A", "Megrivirus A", "Myna chaphamaparvovirus"),
  lib_id = c("SRR11775227", "SRR16071453", "SRR16071454"),
  stringsAsFactors = FALSE
)
false_positive_cells$abundance_after_filter <- vapply(seq_len(nrow(false_positive_cells)), function(i) {
  species <- false_positive_cells$species_name[i]
  library <- false_positive_cells$lib_id[i]
  if (!species %in% rownames(abundance)) stop("Missing false-positive taxon: ", species)
  if (!library %in% colnames(abundance)) stop("Missing false-positive library: ", library)
  abundance[species, library]
}, numeric(1))
if (any(false_positive_cells$abundance_after_filter != 0)) {
  stop("The three confirmed false-positive cells are not all zero.")
}
write.csv(
  false_positive_cells,
  file.path(script_dir, "FigureS3_false_positive_qc.csv"),
  row.names = FALSE
)

annotation <- annotation[match(library_columns, annotation$lib_id), , drop = FALSE]
sample_counts <- table(annotation[[sample_type_col]])

generate_subsets <- function(n_libraries, depth, n_iter = N_ITER) {
  total_combinations <- choose(n_libraries, depth)
  if (is.finite(total_combinations) && total_combinations <= n_iter) {
    return(t(combn(seq_len(n_libraries), depth)))
  }

  draws <- matrix(NA_integer_, nrow = n_iter, ncol = depth)
  for (i in seq_len(n_iter)) {
    draws[i, ] <- sample.int(n_libraries, size = depth, replace = FALSE)
  }
  draws
}

subset_metrics <- function(taxon_by_library, subsets, chunk_size = 500L) {
  richness <- numeric(nrow(subsets))
  shannon <- numeric(nrow(subsets))

  starts <- seq(1L, nrow(subsets), by = chunk_size)
  for (start in starts) {
    stop_i <- min(start + chunk_size - 1L, nrow(subsets))
    selected <- subsets[start:stop_i, , drop = FALSE]
    n_draws <- nrow(selected)
    selector <- matrix(0, nrow = ncol(taxon_by_library), ncol = n_draws)
    for (j in seq_len(ncol(selected))) {
      selector[cbind(selected[, j], seq_len(n_draws))] <- 1
    }

    pooled <- taxon_by_library %*% selector
    positive <- pooled > 0
    richness[start:stop_i] <- colSums(positive)
    totals <- colSums(pooled)
    proportions <- matrix(0, nrow = nrow(pooled), ncol = ncol(pooled))
    nonzero_totals <- totals > 0
    if (any(nonzero_totals)) {
      proportions[, nonzero_totals] <- sweep(
        pooled[, nonzero_totals, drop = FALSE],
        2L, totals[nonzero_totals], "/"
      )
    }
    log_proportions <- matrix(0, nrow = nrow(proportions), ncol = ncol(proportions))
    positive_p <- proportions > 0
    log_proportions[positive_p] <- log(proportions[positive_p])
    shannon[start:stop_i] <- -colSums(proportions * log_proportions)
  }
  list(richness = as.numeric(richness), shannon = as.numeric(shannon))
}

summarize_values <- function(x) {
  data.frame(
    mean = mean(x),
    median = median(x),
    lower_95 = as.numeric(stats::quantile(x, 0.025, names = FALSE)),
    upper_95 = as.numeric(stats::quantile(x, 0.975, names = FALSE)),
    min = min(x),
    max = max(x),
    stringsAsFactors = FALSE
  )
}

set.seed(SEED)
eligible_types <- TYPE_ORDER[
  vapply(TYPE_ORDER, function(x) {
    !is.na(sample_counts[x]) && sample_counts[x] >= TYPE_DEPTH
  }, logical(1))
]

overview_rows <- list()
curve_rows <- list()
type_summary_rows <- list()
n5_draw_rows <- list()

for (sample_type in names(sample_counts)) {
  ids <- annotation$lib_id[annotation[[sample_type_col]] == sample_type]
  observed <- abundance[, ids, drop = FALSE]
  pooled <- rowSums(observed)
  total <- sum(pooled)
  p <- if (total > 0) pooled / total else rep(0, length(pooled))
  positive <- p > 0
  observed_shannon <- if (any(positive)) -sum(p[positive] * log(p[positive])) else 0

  overview_rows[[length(overview_rows) + 1L]] <- data.frame(
    sample_type = sample_type,
    anatomical_class = if (sample_type %in% names(TYPE_CLASS)) TYPE_CLASS[[sample_type]] else "Unclassified",
    n_libraries = length(ids),
    observed_cumulative_richness = sum(pooled > 0),
    observed_pooled_shannon = observed_shannon,
    included_in_n5_analysis = sample_type %in% eligible_types,
    stringsAsFactors = FALSE
  )
}

for (sample_type in eligible_types) {
  ids <- annotation$lib_id[annotation[[sample_type_col]] == sample_type]
  x <- abundance[, ids, drop = FALSE]
  n_libraries <- length(ids)

  for (depth in seq_len(TYPE_DEPTH)) {
    subsets <- generate_subsets(n_libraries, depth)
    metrics <- subset_metrics(x, subsets)

    for (metric_name in c("Taxonomic richness", "Shannon diversity")) {
      values <- if (metric_name == "Taxonomic richness") metrics$richness else metrics$shannon
      stats <- summarize_values(values)
      curve_rows[[length(curve_rows) + 1L]] <- data.frame(
        sample_type = sample_type,
        anatomical_class = TYPE_CLASS[[sample_type]],
        available_libraries = n_libraries,
        rarefaction_depth = depth,
        metric = metric_name,
        n_subsets = length(values),
        stats,
        stringsAsFactors = FALSE
      )
    }

    if (depth == TYPE_DEPTH) {
      n5_draw_rows[[length(n5_draw_rows) + 1L]] <- data.frame(
        sample_type = sample_type,
        anatomical_class = TYPE_CLASS[[sample_type]],
        draw = seq_along(metrics$richness),
        richness = metrics$richness,
        shannon = metrics$shannon,
        stringsAsFactors = FALSE
      )

      richness_stats <- summarize_values(metrics$richness)
      shannon_stats <- summarize_values(metrics$shannon)
      type_summary_rows[[length(type_summary_rows) + 1L]] <- data.frame(
        sample_type = sample_type,
        anatomical_class = TYPE_CLASS[[sample_type]],
        available_libraries = n_libraries,
        rarefaction_depth = TYPE_DEPTH,
        n_subsets = length(metrics$richness),
        richness_mean = richness_stats$mean,
        richness_median = richness_stats$median,
        richness_lower_95 = richness_stats$lower_95,
        richness_upper_95 = richness_stats$upper_95,
        shannon_mean = shannon_stats$mean,
        shannon_median = shannon_stats$median,
        shannon_lower_95 = shannon_stats$lower_95,
        shannon_upper_95 = shannon_stats$upper_95,
        stringsAsFactors = FALSE
      )
    }
  }
}

overview <- do.call(rbind, overview_rows)
type_summary <- do.call(rbind, type_summary_rows)
type_summary$plot_order <- match(type_summary$sample_type, TYPE_ORDER)
type_summary <- type_summary[order(type_summary$plot_order), , drop = FALSE]
type_summary$plot_order <- NULL
curves <- do.call(rbind, curve_rows)
n5_draws <- do.call(rbind, n5_draw_rows)

category_ids <- lapply(names(CLASS_COLORS)[1:2], function(category) {
  types <- names(TYPE_CLASS)[TYPE_CLASS == category]
  annotation$lib_id[annotation[[sample_type_col]] %in% types]
})
names(category_ids) <- names(CLASS_COLORS)[1:2]
category_depth <- min(vapply(category_ids, length, integer(1)))

category_rows <- list()
category_draw_rows <- list()
for (category in names(category_ids)) {
  ids <- category_ids[[category]]
  x <- abundance[, ids, drop = FALSE]

  for (depth in seq_len(category_depth)) {
    subsets <- generate_subsets(length(ids), depth)
    metrics <- subset_metrics(x, subsets)

    for (metric_name in c("Taxonomic richness", "Shannon diversity")) {
      values <- if (metric_name == "Taxonomic richness") metrics$richness else metrics$shannon
      stats <- summarize_values(values)
      category_rows[[length(category_rows) + 1L]] <- data.frame(
        anatomical_class = category,
        available_libraries = length(ids),
        rarefaction_depth = depth,
        metric = metric_name,
        n_subsets = length(values),
        stats,
        stringsAsFactors = FALSE
      )
    }

    if (depth == category_depth) {
      category_draw_rows[[length(category_draw_rows) + 1L]] <- data.frame(
        anatomical_class = category,
        draw = seq_along(metrics$richness),
        richness = metrics$richness,
        shannon = metrics$shannon,
        stringsAsFactors = FALSE
      )
    }
  }
}

category_curves <- do.call(rbind, category_rows)
category_draws <- do.call(rbind, category_draw_rows)

write.csv(overview, file.path(script_dir, "FigureS3_sample_type_overview_R.csv"), row.names = FALSE)
write.csv(type_summary, file.path(script_dir, "FigureS3_rarefaction_n5_summary_R.csv"), row.names = FALSE)
write.csv(curves, file.path(script_dir, "FigureS3_rarefaction_curve_summary_by_sample_type_R.csv"), row.names = FALSE)
write.csv(n5_draws, file.path(script_dir, "FigureS3_rarefaction_draws_n5_R.csv"), row.names = FALSE)
write.csv(category_curves, file.path(script_dir, "FigureS3_rarefaction_curve_summary_by_anatomical_class_R.csv"), row.names = FALSE)
write.csv(category_draws, file.path(script_dir, "FigureS3_rarefaction_draws_matched_anatomical_class_R.csv"), row.names = FALSE)

matched_summary <- category_curves[
  category_curves$metric == "Taxonomic richness" &
    category_curves$rarefaction_depth == category_depth, ,
  drop = FALSE
]
write.csv(matched_summary, file.path(script_dir, "FigureS3_matched_anatomical_class_summary_R.csv"), row.names = FALSE)

grid_col <- "#D9D9D9"
point_colors <- unname(CLASS_COLORS[type_summary$anatomical_class])

plot_panel_A <- function() {
  d <- category_curves[category_curves$metric == "Taxonomic richness", , drop = FALSE]
  y_top <- ceiling(max(d$upper_95, na.rm = TRUE) / 25) * 25
  x_ticks <- seq(0, category_depth, by = 5)
  y_ticks <- seq(0, y_top, by = 25)

  plot(NA, NA, xlim = c(0, category_depth), ylim = c(0, max(y_ticks)),
       axes = FALSE, xlab = "", ylab = "", bty = "n")
  abline(h = y_ticks, col = grid_col, lwd = 1)

  for (category in names(CLASS_COLORS)[1:2]) {
    z <- d[d$anatomical_class == category, , drop = FALSE]
    z <- z[order(z$rarefaction_depth), , drop = FALSE]
    origin <- z[1, , drop = FALSE]
    origin$rarefaction_depth <- 0
    origin$mean <- 0
    origin$lower_95 <- 0
    origin$upper_95 <- 0
    z <- rbind(origin, z)
    polygon(
      c(z$rarefaction_depth, rev(z$rarefaction_depth)),
      c(z$lower_95, rev(z$upper_95)),
      col = grDevices::adjustcolor(CLASS_COLORS[[category]], alpha.f = 0.16),
      border = NA
    )
    lines(z$rarefaction_depth, z$mean, col = CLASS_COLORS[[category]], lwd = 1.8)
  }

  axis(1, at = x_ticks, labels = x_ticks)
  axis(2, at = y_ticks, labels = y_ticks, las = 1)
  mtext("Number of libraries", side = 1, line = 2.4)
  mtext("Expected taxonomic richness", side = 2, line = 2.8)

  legend(
    "topleft", legend = names(CLASS_COLORS)[1:2],
    col = unname(CLASS_COLORS[1:2]), lwd = 1.8,
    bty = "n", horiz = TRUE, cex = 0.78,
    seg.len = 2.2, x.intersp = 0.7
  )
  mtext("A", side = 3, adj = 0, line = 0.1, cex = 1.15, font = 2)
}

plot_forest_panel <- function(panel_id, median_col, lower_col, upper_col,
                              xlab, letter, show_y_axis = TRUE) {
  x_values <- type_summary[[median_col]]
  lower <- type_summary[[lower_col]]
  upper <- type_summary[[upper_col]]
  y <- rev(seq_len(nrow(type_summary))) - 1L

  x_max <- max(upper, na.rm = TRUE)
  x_ticks <- if (panel_id == "B") seq(0, 120, by = 20) else seq(0, 3.0, by = 0.5)
  x_lim <- c(0, max(x_ticks, x_max * 1.08))

  plot(NA, NA, xlim = x_lim, ylim = c(-0.20, nrow(type_summary) - 0.35),
       axes = FALSE, xlab = "", ylab = "", bty = "n")
  abline(v = x_ticks, col = grid_col, lwd = 1)
  axis(1, at = x_ticks, labels = x_ticks)
  if (show_y_axis) {
    axis(2, at = y, labels = type_summary$sample_type, las = 1, cex.axis = 0.76)
  } else {
    axis(2, at = y, labels = FALSE, tick = FALSE)
  }
  mtext(xlab, side = 1, line = 1.6)

  for (i in seq_along(y)) {
    segments(lower[i], y[i], upper[i], y[i],
             col = point_colors[i], lwd = 1.35, lend = 1)
    points(x_values[i], y[i], pch = 21, bg = point_colors[i],
           col = "white", lwd = 0.55, cex = 1.15)
  }
  mtext(letter, side = 3, adj = 0, line = 0.1, cex = 1.15, font = 2)
}

draw_figure <- function() {
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE),
         heights = c(1.0, 1.2), widths = c(1.08, 0.92))
  par(family = "sans", xpd = FALSE, oma = c(0, 0, 0, 0))

  par(mar = c(3.9, 4.8, 1.1, 0.8))
  plot_panel_A()

  par(mar = c(3.1, 6.3, 1.1, 0.0))
  plot_forest_panel(
    "B", "richness_median", "richness_lower_95", "richness_upper_95",
    "Taxonomic richness", "B", show_y_axis = TRUE
  )

  par(mar = c(3.1, 0.8, 1.1, 0.2))
  plot_forest_panel(
    "C", "shannon_median", "shannon_lower_95", "shannon_upper_95",
    "Shannon diversity", "C", show_y_axis = FALSE
  )

}

pdf_file <- file.path(script_dir, "Figure_S3_rarefaction_R.pdf")
grDevices::pdf(pdf_file, width = 7.2, height = 6.7, useDingbats = FALSE)
draw_figure()
grDevices::dev.off()

png_file <- file.path(script_dir, "Figure_S3_rarefaction_R.png")
grDevices::png(png_file, width = 2160, height = 2010, res = 300, type = "cairo")
draw_figure()
grDevices::dev.off()

tif_file <- file.path(script_dir, "Figure_S3_rarefaction_R.tif")
grDevices::tiff(tif_file, width = 4320, height = 4020, res = 600, compression = "lzw")
draw_figure()
grDevices::dev.off()


