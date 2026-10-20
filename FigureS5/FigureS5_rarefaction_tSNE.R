#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2", "readxl", "Rtsne", "vegan")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(Rtsne)
  library(vegan)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

raref_matrix_file <- file.path(script_dir, "FigureS5_virus_table_filtered.csv")
raref_annotation_file <- file.path(script_dir, "FigureS5_libs_annot.xlsx")
tsne_matrix_file <- file.path(script_dir, "FigureS5_virus_table_selected_samples_updated.csv")
tsne_annotation_file <- file.path(script_dir, "FigureS5_library_annotation.csv")

required_inputs <- c(raref_matrix_file, raref_annotation_file,
                     tsne_matrix_file, tsne_annotation_file)
if (any(!file.exists(required_inputs))) {
  stop("Missing FigureS5 input file(s): ",
       paste(basename(required_inputs[!file.exists(required_inputs)]), collapse = ", "))
}

clean_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\uFEFF", "", x, fixed = TRUE)
  x <- gsub("[\r\n\t]", "", x)
  trimws(x)
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  formatC(p, format = "f", digits = 3)
}

group_colors <- c(
  Charadriiformes = "#1B9E77",
  Passeriformes = "#D95F02",
  Psittaciformes = "#7570B3",
  Anseriformes = "#E7298A",
  Galliformes = "#66A61E",
  Procellariiformes = "#E6AB02",
  China = "#1F78B4",
  `New Zealand` = "#6A3D9A",
  Australia = "#FF7F00",
  `the United States` = "#A6CEE3",
  Wild = "#E31A1C",
  Domestic = "#B15928"
)

target_sample_type <- "Fecal samples"
seed <- 20260902L
curve_replicates <- 1000L
minimum_libraries <- 5L
tsne_seed <- 123L
tsne_perplexity <- 7L
tsne_max_iter <- 3000L
permutations <- 999L

analysis_specs <- list(
  list(id = "Wild-domestic", variable = "Domestic or Wild",
       candidates = c("Wild", "Domestic")),
  list(id = "Host order", variable = "Order",
       candidates = c("Charadriiformes", "Passeriformes", "Psittaciformes",
                      "Anseriformes", "Galliformes", "Procellariiformes")),
  list(id = "Geography", variable = "Location",
       candidates = c("China", "New Zealand", "Australia", "the United States"))
)

## ---------- Rarefaction data and curves ----------
abundance <- read.csv(raref_matrix_file, check.names = FALSE,
                      stringsAsFactors = FALSE, na.strings = c("", "NA"))
annotation <- as.data.frame(
  readxl::read_excel(raref_annotation_file, sheet = 1,
                     .name_repair = "minimal"),
  stringsAsFactors = FALSE
)

if (!"species_name" %in% names(abundance)) {
  stop("The rarefaction matrix must contain a species_name column.")
}
required_annotation <- c("lib_id", "Location", "Domestic or Wild",
                          "Order", "Sample Type")
missing_annotation <- setdiff(required_annotation, names(annotation))
if (length(missing_annotation) > 0L) {
  stop("Missing annotation columns: ", paste(missing_annotation, collapse = ", "))
}

feature_start <- match("species_name", names(abundance)) + 1L
raref_library_ids <- clean_id(names(abundance)[feature_start:ncol(abundance)])
names(abundance)[feature_start:ncol(abundance)] <- raref_library_ids
abundance$species_name <- clean_id(abundance$species_name)
annotation$lib_id <- clean_id(annotation$lib_id)

if (anyDuplicated(raref_library_ids) || anyDuplicated(annotation$lib_id)) {
  stop("Duplicate library IDs found in the rarefaction inputs.")
}
numeric_abundance <- as.matrix(data.frame(lapply(
  abundance[raref_library_ids], function(x) suppressWarnings(as.numeric(x))
), check.names = FALSE))
if (anyDuplicated(abundance$species_name)) {
  if (anyNA(numeric_abundance) || any(numeric_abundance < 0)) {
    stop("Rarefaction abundance values must be non-negative numeric values.")
  }
  numeric_abundance <- rowsum(numeric_abundance, abundance$species_name,
                              reorder = FALSE)
} else {
  rownames(numeric_abundance) <- abundance$species_name
}
if (anyNA(numeric_abundance) || any(numeric_abundance < 0)) {
  stop("Rarefaction abundance values must be non-negative numeric values.")
}

missing_in_annotation <- setdiff(raref_library_ids, annotation$lib_id)
if (length(missing_in_annotation) > 0L) {
  stop("Libraries missing from rarefaction annotation: ",
       paste(head(missing_in_annotation, 10), collapse = ", "))
}
annotation <- annotation[match(raref_library_ids, annotation$lib_id), , drop = FALSE]
annotation$Sample_Type_Clean <- trimws(as.character(annotation$`Sample Type`))
presence_matrix <- numeric_abundance > 0

get_group_info <- function(spec) {
  values <- trimws(as.character(annotation[[spec$variable]]))
  keep <- annotation$Sample_Type_Clean == target_sample_type
  counts <- table(values[keep])
  eligible <- spec$candidates[
    spec$candidates %in% names(counts) &
      as.integer(counts[spec$candidates]) >= minimum_libraries
  ]
  if (length(eligible) < 2L) {
    stop("Fewer than two eligible groups for ", spec$id, ".")
  }
  indices <- lapply(eligible, function(g) which(values == g & keep))
  names(indices) <- eligible
  list(variable = values, groups = eligible, indices = indices, counts = counts)
}

sample_group_presence <- function(indices, groups, depth) {
  sampled <- lapply(indices[groups], function(x) sample(x, depth, replace = FALSE))
  out <- do.call(cbind, lapply(sampled, function(x) {
    rowSums(presence_matrix[, x, drop = FALSE]) > 0
  }))
  colnames(out) <- groups
  out
}

simulate_curve <- function(spec, info) {
  max_depth <- min(vapply(info$indices[info$groups], length, integer(1)))
  result <- vector("list", max_depth * curve_replicates)
  z <- 1L
  for (depth in seq_len(max_depth)) {
    for (replicate_id in seq_len(curve_replicates)) {
      present <- sample_group_presence(info$indices, info$groups, depth)
      richness <- colSums(present)
      result[[z]] <- data.frame(
        analysis = spec$id, group = names(richness), depth = depth,
        replicate = replicate_id, richness = as.numeric(richness),
        stringsAsFactors = FALSE
      )
      z <- z + 1L
    }
  }
  do.call(rbind, result)
}

set.seed(seed)
curve_raw <- vector("list", length(analysis_specs))
excluded_rows <- vector("list", length(analysis_specs))
settings_rows <- vector("list", length(analysis_specs))
group_info_by_analysis <- vector("list", length(analysis_specs))
names(group_info_by_analysis) <- vapply(analysis_specs, `[[`, character(1), "id")

for (i in seq_along(analysis_specs)) {
  spec <- analysis_specs[[i]]
  info <- get_group_info(spec)
  group_info_by_analysis[[spec$id]] <- info
  curve_raw[[i]] <- simulate_curve(spec, info)
  counts <- as.integer(info$counts)
  names(counts) <- names(info$counts)
  excluded <- spec$candidates[
    !(spec$candidates %in% names(counts)) |
      as.integer(counts[spec$candidates]) < minimum_libraries
  ]
  if (spec$id == "Host order") excluded <- unique(c(excluded, "pooling"))
  excluded_rows[[i]] <- data.frame(
    analysis = rep(spec$id, length(excluded)), group = excluded,
    library_count = unname(counts[excluded]),
    reason = ifelse(spec$id == "Host order" & excluded == "pooling",
                    "not a biological host-order stratum",
                    paste0("fewer than ", minimum_libraries,
                           " libraries in ", target_sample_type)),
    stringsAsFactors = FALSE
  )
  settings_rows[[i]] <- data.frame(
    analysis = spec$id,
    eligible_groups = paste(info$groups, collapse = "; "),
    sample_type = target_sample_type,
    libraries_per_group = min(vapply(info$indices[info$groups], length, integer(1))),
    curve_replicates = curve_replicates, presence_rule = "abundance > 0",
    stringsAsFactors = FALSE
  )
}

curve_raw <- do.call(rbind, curve_raw)
curve_key <- interaction(curve_raw$analysis, curve_raw$group,
                          curve_raw$depth, drop = TRUE)
curve_summary <- do.call(rbind, lapply(split(curve_raw, curve_key), function(x) {
  data.frame(
    analysis = x$analysis[1], group = x$group[1], depth = x$depth[1],
    mean_richness = mean(x$richness),
    lower_95 = unname(quantile(x$richness, 0.025)),
    upper_95 = unname(quantile(x$richness, 0.975)),
    stringsAsFactors = FALSE
  )
}))

analysis_levels <- vapply(analysis_specs, `[[`, character(1), "id")
curve_summary$analysis <- factor(curve_summary$analysis, levels = analysis_levels)

make_rarefaction_plot <- function(analysis_id, panel_label) {
  dat <- curve_summary[as.character(curve_summary$analysis) == analysis_id, , drop = FALSE]
  groups <- group_info_by_analysis[[analysis_id]]$groups
  dat$group <- factor(dat$group, levels = groups)
  max_depth <- max(dat$depth)
  labels <- do.call(rbind, lapply(split(dat, dat$group), function(x) {
    x[which.max(x$depth), , drop = FALSE]
  }))
  labels$label_x <- labels$depth + max(0.12, max_depth * 0.025)
  ggplot(dat, aes(depth, mean_richness, colour = group, fill = group,
                  group = group)) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.65) +
    geom_text(data = labels,
              aes(x = label_x, y = mean_richness, label = group, colour = group),
              inherit.aes = FALSE, hjust = 0, size = 2.15, show.legend = FALSE) +
    annotate("text", x = -Inf, y = Inf, label = panel_label,
             hjust = -0.18, vjust = 1.2, size = 4.8, fontface = "bold") +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_x_continuous(breaks = seq(0, max_depth, by = 2),
                       limits = c(0, max_depth + max(0.4, max_depth * 0.08)),
                       expand = c(0, 0)) +
    labs(x = "Libraries per group", y = "Detected taxa per group") +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 8.5, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.text = element_text(colour = "black", size = 6.5),
      axis.title = element_text(colour = "black", size = 7.5),
      legend.position = "none", plot.margin = margin(7, 25, 5, 8)
    )
}

rarefaction_plots <- list(
  make_rarefaction_plot("Wild-domestic", "A"),
  make_rarefaction_plot("Host order", "B"),
  make_rarefaction_plot("Geography", "C")
)

## ---------- tSNE plots matching the three rarefaction columns ----------
tsne_abundance <- read.csv(tsne_matrix_file, check.names = FALSE,
                           stringsAsFactors = FALSE, na.strings = c("", "NA"))
tsne_annotation <- read.csv(tsne_annotation_file, check.names = FALSE,
                            stringsAsFactors = FALSE, na.strings = c("", "NA"))
if (!"species_name" %in% names(tsne_abundance)) {
  stop("The tSNE matrix must contain a species_name column.")
}
required_tsne_annotation <- c("lib_id", "Location", "Domestic or Wild",
                              "Order", "Sample Type")
if (length(setdiff(required_tsne_annotation, names(tsne_annotation))) > 0L) {
  stop("The tSNE annotation is missing required columns.")
}

tsne_library_ids <- clean_id(setdiff(names(tsne_abundance), "species_name"))
names(tsne_abundance)[match(tsne_library_ids, names(tsne_abundance))] <- tsne_library_ids
tsne_abundance$species_name <- clean_id(tsne_abundance$species_name)
tsne_annotation$lib_id <- clean_id(tsne_annotation$lib_id)
tsne_matrix <- as.matrix(data.frame(lapply(
  tsne_abundance[tsne_library_ids], function(x) suppressWarnings(as.numeric(x))
), check.names = FALSE))
if (anyNA(tsne_matrix) || any(tsne_matrix < 0)) {
  stop("tSNE abundance values must be non-negative numeric values.")
}
if (anyDuplicated(tsne_abundance$species_name)) {
  tsne_matrix <- rowsum(tsne_matrix, tsne_abundance$species_name, reorder = FALSE)
} else {
  rownames(tsne_matrix) <- tsne_abundance$species_name
}
colnames(tsne_matrix) <- tsne_library_ids

run_tsne_analysis <- function(analysis_id, panel_label) {
  info <- group_info_by_analysis[[analysis_id]]
  spec <- analysis_specs[[match(analysis_id, analysis_levels)]]
  variable <- trimws(as.character(tsne_annotation[[spec$variable]]))
  sample_type <- trimws(as.character(tsne_annotation$`Sample Type`))
  keep_meta <- sample_type == target_sample_type & variable %in% info$groups
  target_ids <- intersect(tsne_annotation$lib_id[keep_meta], colnames(tsne_matrix))
  if (length(target_ids) < 4L) stop("Too few tSNE libraries for ", analysis_id, ".")
  x <- tsne_matrix[, target_ids, drop = FALSE]
  x <- x[, colSums(x) > 0, drop = FALSE]
  if (ncol(x) < 4L) stop("Too few non-zero tSNE libraries for ", analysis_id, ".")
  meta <- tsne_annotation[match(colnames(x), tsne_annotation$lib_id), , drop = FALSE]
  meta$group <- factor(trimws(as.character(meta[[spec$variable]])), levels = info$groups)
  if (anyNA(meta$group) || nlevels(droplevels(meta$group)) < 2L) {
    stop("At least two groups are required for tSNE analysis: ", analysis_id)
  }
  x <- log10(x + 1)
  bray <- vegan::vegdist(t(x), method = "bray")
  ad <- vegan::adonis2(bray ~ group, data = meta, permutations = permutations,
                       by = "terms")
  bd <- vegan::betadisper(bray, meta$group, add = "lingoes")
  bd_test <- vegan::permutest(bd, permutations = permutations)
  r2 <- unname(ad$R2[1]); perm_p <- unname(ad$`Pr(>F)`[1])
  disp_p <- unname(bd_test$tab[1, "Pr(>F)"])
  n <- nrow(meta)
  perplexity <- min(tsne_perplexity, floor((n - 1) / 3))
  if (perplexity < 2L) stop("Too few libraries for a stable tSNE perplexity.")
  set.seed(tsne_seed)
  fit <- Rtsne::Rtsne(as.matrix(bray), is_distance = TRUE,
                      perplexity = perplexity, max_iter = tsne_max_iter,
                      check_duplicates = FALSE, pca = FALSE, verbose = FALSE)
  coords <- data.frame(
    analysis = analysis_id, panel = panel_label, lib_id = meta$lib_id,
    group = as.character(meta$group), tSNE1 = fit$Y[, 1], tSNE2 = fit$Y[, 2],
    stringsAsFactors = FALSE
  )
  stats <- data.frame(
    analysis = analysis_id, panel = panel_label, n_libraries = n,
    groups = paste(info$groups, collapse = "; "), perplexity = perplexity,
    seed = tsne_seed, permanova_R2 = r2, permanova_P = perm_p,
    permdisp_P = disp_p, stringsAsFactors = FALSE
  )
  list(coords = coords, stats = stats)
}

tsne_runs <- list(
  run_tsne_analysis("Wild-domestic", "D"),
  run_tsne_analysis("Host order", "E"),
  run_tsne_analysis("Geography", "F")
)
names(tsne_runs) <- analysis_levels

make_tsne_plot <- function(run) {
  dat <- run$coords; st <- run$stats
  x_range <- range(dat$tSNE1); y_range <- range(dat$tSNE2)
  x_span <- max(diff(x_range), 1)
  x_pad <- max(x_span * 0.04, 0.1)
  right_pad <- max(x_span * 0.28, 0.1)
  y_pad <- max(diff(y_range) * 0.05, 0.1)
  stat_x <- x_range[2] + right_pad * 0.90
  stat_label <- paste0(
    "PERMANOVA: R² = ", formatC(st$permanova_R2, format = "f", digits = 3),
    ", P = ", format_p(st$permanova_P), "\n",
    "PERMDISP: P = ", format_p(st$permdisp_P)
  )
  ggplot(dat, aes(tSNE1, tSNE2, colour = group)) +
    geom_point(size = 2.1, alpha = 0.72) +
    annotate("text", x = -Inf, y = Inf, label = st$panel,
             hjust = -0.18, vjust = 1.2, size = 4.8, fontface = "bold") +
    annotate("text", x = stat_x, y = y_range[2] - y_pad,
             label = stat_label, hjust = 1, vjust = 1, size = 2.05,
             colour = "black") +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    labs(x = "t-SNE1", y = "t-SNE2", colour = NULL) +
    # Do not force a 1:1 physical aspect ratio here: the three tSNE panels
    # have different coordinate ranges, and coord_equal() would make their
    # visible plotting widths different inside the shared three-column grid.
    coord_cartesian(
      xlim = c(x_range[1] - x_pad, x_range[2] + right_pad),
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
      expand = FALSE, clip = "off"
    ) +
    theme_minimal(base_size = 8.5, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.25),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.text = element_text(colour = "black", size = 6.3),
      axis.title = element_text(colour = "black", size = 7.5),
      legend.position = "bottom", legend.direction = "horizontal",
      legend.text = element_text(size = 5.7),
      legend.key.height = grid::unit(0.22, "cm"),
      legend.key.width = grid::unit(0.32, "cm"), plot.margin = margin(7, 8, 2, 8)
    )
}

tsne_plots <- lapply(tsne_runs, make_tsne_plot)

## ---------- Machine-readable outputs ----------
write.csv(curve_summary, file.path(script_dir,
                                   "FigureS5_BCD_group_richness_rarefaction.csv"),
          row.names = FALSE, na = "")
write.csv(do.call(rbind, excluded_rows), file.path(script_dir,
                                                   "FigureS5_BCD_excluded_groups.csv"),
          row.names = FALSE, na = "")
write.csv(do.call(rbind, settings_rows), file.path(script_dir,
                                                  "FigureS5_BCD_analysis_settings.csv"),
          row.names = FALSE, na = "")
write.csv(do.call(rbind, lapply(tsne_runs, `[[`, "coords")),
          file.path(script_dir, "FigureS5_tSNE_coordinates.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(tsne_runs, `[[`, "stats")),
          file.path(script_dir, "FigureS5_tSNE_statistics.csv"), row.names = FALSE)

summary_lines <- c(
  "Figure S5: rarefaction curves with matched tSNE panels",
  "Top row: A Wild-domestic, B Host order, C Geography",
  "Bottom row: D Wild-domestic, E Host order, F Geography",
  paste0("Rarefaction sample type: ", target_sample_type),
  paste0("Rarefaction replicates: ", curve_replicates),
  paste0("tSNE transform: log10(abundance + 1); maximum iterations: ", tsne_max_iter),
  paste0("tSNE seed: ", tsne_seed, "; requested perplexity: ", tsne_perplexity),
  paste0("PERMANOVA/PERMDISP permutations: ", permutations),
  "Each tSNE panel uses the eligible groups from its matching rarefaction panel."
)
writeLines(summary_lines, file.path(script_dir, "FigureS5_README.txt"))

## ---------- Combined six-panel export ----------
panel_grobs <- lapply(c(rarefaction_plots, tsne_plots), ggplotGrob)
draw_combined <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(2, 3,
                               heights = grid::unit(c(1, 1.06), "null"),
                               widths = grid::unit(rep(1, 3), "null")),
    gp = grid::gpar(fontfamily = "sans")
  ))
  for (i in seq_along(panel_grobs)) {
    row <- if (i <= 3L) 1L else 2L
    col <- ((i - 1L) %% 3L) + 1L
    grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
    grid::grid.draw(panel_grobs[[i]])
    grid::popViewport()
  }
  grid::popViewport()
}

pdf(file.path(script_dir, "FigureS5_rarefaction_tSNE_combined.pdf"),
    width = 11.2, height = 7.6, useDingbats = FALSE)
draw_combined()
dev.off()

png(file.path(script_dir, "FigureS5_rarefaction_tSNE_combined.png"),
    width = 3360, height = 2280, res = 300, type = "cairo")
draw_combined()
dev.off()

message("FigureS5 combined rarefaction+tSNE figure completed: ", script_dir)
