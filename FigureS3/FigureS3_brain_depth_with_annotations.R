#!/usr/bin/env Rscript

# Figure S3: brain-positive depth profiles with NR-supported protein
# annotation tracks.  The protein tracks use gggenes and share the same
# nucleotide coordinate system as the depth profiles.

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2", "gggenes", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install required R package(s) before running: ",
    paste(missing_packages, collapse = ", "),
    "\nFor example: install.packages(c(\"ggplot2\", \"dplyr\")); ",
    "install.packages(\"gggenes\")"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(gggenes)
  library(dplyr)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  getwd()
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[1], fixed = TRUE)
}

input_dir <- get_arg("input-dir", script_dir)
outdir <- get_arg("outdir", script_dir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(input_dir, "FigureS3_brain_depth_summary.tsv")
normalization_file <- file.path(input_dir, "FigureS3_library_normalization.tsv")
targets_file <- file.path(input_dir, "FigureS3_four_reference_targets.tsv")
annotation_file <- file.path(
  input_dir,
  "FigureS3_NR_NCBI_annotation_segments_for_plot.tsv"
)

for (required_file in c(
  summary_file,
  normalization_file,
  targets_file,
  annotation_file
)) {
  if (!file.exists(required_file)) {
    stop("Required input file not found: ", required_file)
  }
}

depth_dir <- file.path(input_dir, "depth")
if (!dir.exists(depth_dir)) {
  # The downloaded files may be placed directly in FigureS3.
  depth_dir <- input_dir
}

summary_data <- read.delim(
  summary_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
if ("status" %in% names(summary_data)) {
  summary_data <- summary_data[summary_data$status == "completed", , drop = FALSE]
}
if (nrow(summary_data) == 0) {
  stop("No completed depth profiles were found in the summary table.")
}

required_summary_columns <- c(
  "species_name", "library", "rpm", "breadth_1x", "bam_reference"
)
missing_summary_columns <- setdiff(required_summary_columns, names(summary_data))
if (length(missing_summary_columns) > 0) {
  stop(
    "Summary table is missing column(s): ",
    paste(missing_summary_columns, collapse = ", ")
  )
}

normalization_data <- read.delim(
  normalization_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_normalization_columns <- c("library", "norRNA_reads")
missing_normalization_columns <- setdiff(
  required_normalization_columns,
  names(normalization_data)
)
if (length(missing_normalization_columns) > 0) {
  stop(
    "Normalization table is missing column(s): ",
    paste(missing_normalization_columns, collapse = ", ")
  )
}

normalization_data$library <- as.character(normalization_data$library)
normalization_data$norRNA_reads <- as.numeric(normalization_data$norRNA_reads)
summary_data$library <- as.character(summary_data$library)
summary_data$norRNA_reads <- normalization_data$norRNA_reads[
  match(summary_data$library, normalization_data$library)
]
if (anyNA(summary_data$norRNA_reads)) {
  missing_libraries <- unique(summary_data$library[is.na(summary_data$norRNA_reads)])
  stop(
    "Missing norRNA_reads for library/libraries: ",
    paste(missing_libraries, collapse = ", ")
  )
}
summary_data$reads <- round(
  as.numeric(summary_data$rpm) * summary_data$norRNA_reads / 1e6
)

targets <- read.delim(
  targets_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_target_columns <- c("sequence_id", "bam_reference")
missing_target_columns <- setdiff(required_target_columns, names(targets))
if (length(missing_target_columns) > 0) {
  stop(
    "Target table is missing column(s): ",
    paste(missing_target_columns, collapse = ", ")
  )
}

annotation_data <- read.delim(
  annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
required_annotation_columns <- c(
  "sequence_id", "strand", "retained_start", "retained_end",
  "putative_gene_product", "accession", "identity_pct", "evalue"
)
missing_annotation_columns <- setdiff(
  required_annotation_columns,
  names(annotation_data)
)
if (length(missing_annotation_columns) > 0) {
  stop(
    "Annotation table is missing column(s): ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

# Keep all annotation coordinates in the retained reference coordinate system.
annotation_data <- annotation_data %>%
  mutate(
    start = as.numeric(retained_start),
    end = as.numeric(retained_end),
    forward = as.character(strand) == "+",
    putative_gene_product = as.character(putative_gene_product),
    accession = as.character(accession),
    annotation_class = case_when(
      grepl("hexon|capsid|pVI", putative_gene_product, ignore.case = TRUE) ~
        "Structural protein",
      grepl("replicase|polymerase|protease", putative_gene_product, ignore.case = TRUE) ~
        "Replication-associated protein",
      TRUE ~ "Other viral protein"
    ),
    label = case_when(
      grepl("hexon", putative_gene_product, ignore.case = TRUE) ~ "hexon",
      grepl("pVI", putative_gene_product, ignore.case = TRUE) ~ "pVI",
      grepl("capsid", putative_gene_product, ignore.case = TRUE) ~ "capsid",
      grepl("replicase", putative_gene_product, ignore.case = TRUE) ~ "replicase",
      grepl("large polymerase", putative_gene_product, ignore.case = TRUE) ~ "L polymerase",
      grepl("polymerase", putative_gene_product, ignore.case = TRUE) ~ "polymerase",
      grepl("protease", putative_gene_product, ignore.case = TRUE) ~ "protease",
      TRUE ~ putative_gene_product
    )
  ) %>%
  filter(is.finite(start), is.finite(end), end >= start)

depth_files <- list.files(
  depth_dir,
  pattern = "\\.depth\\.tsv$",
  full.names = TRUE
)
if (length(depth_files) == 0) {
  stop("No *.depth.tsv files found in: ", depth_dir)
}

normalise_path <- function(x) gsub("\\\\", "/", x)

find_depth_file <- function(summary_row, index) {
  if ("depth_file" %in% names(summary_row)) {
    requested <- basename(normalise_path(as.character(summary_row$depth_file)))
    exact <- depth_files[basename(depth_files) == requested]
    if (length(exact) == 1) return(exact)
  }
  if (index <= length(depth_files)) return(depth_files[index])
  NA_character_
}

wrap_title <- function(x, width = 42) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

pretty_depth_scale <- function(max_depth) {
  if (!is.finite(max_depth) || max_depth <= 0) {
    return(list(breaks = c(0, 1), upper = 1))
  }
  breaks <- pretty(c(0, max_depth), n = 5)
  breaks <- breaks[is.finite(breaks) & breaks >= 0]
  upper <- max(breaks[breaks >= max_depth], na.rm = TRUE)
  if (!is.finite(upper) || upper <= max_depth) {
    step <- if (length(breaks) >= 2) diff(tail(breaks, 2)) else max(1, max_depth / 5)
    upper <- max(max(breaks), max_depth + step)
    breaks <- sort(unique(c(breaks, upper)))
  }
  list(breaks = breaks[breaks <= upper], upper = upper)
}

annotation_colors <- c(
  "Structural protein" = "#4C78A8",
  "Replication-associated protein" = "#E45756",
  "Other viral protein" = "#72B7B2"
)

make_depth_annotation_plot <- function(
    depth_data,
    species_name,
    library_id,
    reads,
    rpm,
    coverage_pct,
    panel_label,
    sequence_id,
    panel_annotations
) {
  max_depth <- max(depth_data$depth, na.rm = TRUE)
  depth_scale <- pretty_depth_scale(max_depth)
  y_upper <- depth_scale$upper

  # Reserve a proportional lower strip for the protein track. The depth axis
  # remains non-negative; the negative range is used only for annotations.
  gene_track_height <- max(y_upper * 0.14, 2)
  gene_y <- -gene_track_height * 0.55

  panel_annotations <- panel_annotations %>%
    mutate(
      gene_y = gene_y
    )

  plot_title <- wrap_title(paste0(species_name, " | ", library_id))
  stats_label <- paste0(
    "reads = ", formatC(as.integer(round(as.numeric(reads))), format = "d", big.mark = ","),
    "\nRPM = ", formatC(as.numeric(rpm), format = "f", digits = 3),
    "\nCoverage = ", formatC(as.numeric(coverage_pct), format = "f", digits = 1), "%"
  )

  ggplot(depth_data, aes(x = position, y = depth)) +
    geom_area(fill = "#648AB4", alpha = 1, colour = NA) +
    geom_line(colour = "#315B88", linewidth = 0.35) +
    geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.35) +
    gggenes::geom_gene_arrow(
      data = panel_annotations,
      aes(
        xmin = start,
        xmax = end,
        y = gene_y,
        forward = forward,
        fill = annotation_class
      ),
      inherit.aes = FALSE,
      arrowhead_height = grid::unit(3, "mm"),
      arrowhead_width = grid::unit(2, "mm"),
      arrow_body_height = grid::unit(2, "mm"),
      colour = "white",
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    # Keep labels horizontal and centred inside their corresponding gene
    # arrows, rather than placing them above the track.
    geom_text(
      data = panel_annotations,
      aes(x = (start + end) / 2, y = gene_y, label = label),
      inherit.aes = FALSE,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 2.35,
      fontface = "bold",
      colour = "#222222",
      check_overlap = TRUE
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = stats_label,
      hjust = 1.05,
      vjust = 1.15,
      size = 2.65
    ) +
    labs(
      title = plot_title,
      tag = panel_label,
      x = "Position (bp)",
      y = "Depth"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(
      breaks = depth_scale$breaks,
      labels = depth_scale$breaks,
      limits = c(-gene_track_height, y_upper),
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = annotation_colors, drop = FALSE) +
    coord_cartesian(
      ylim = c(-gene_track_height, y_upper),
      expand = FALSE,
      clip = "off"
    ) +
    # Use the system sans-serif family. This avoids Windows font-database
    # warnings when Arial is not registered in the current R session.
    theme_classic(base_family = "sans", base_size = 8.5) +
    theme(
      plot.title = element_text(
        size = 8.5,
        hjust = 0.5,
        lineheight = 0.95,
        face = "italic"
      ),
      plot.tag = element_text(size = 11, face = "bold", hjust = 0, vjust = 1),
      plot.tag.position = c(0.01, 0.99),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.2, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      panel.grid.major.y = element_line(
        colour = "grey85", linewidth = 0.3, linetype = "dashed"
      ),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(5, 7, 8, 8)
    )
}

plot_list <- vector("list", nrow(summary_data))
for (i in seq_len(nrow(summary_data))) {
  depth_file <- find_depth_file(summary_data[i, , drop = FALSE], i)
  if (is.na(depth_file) || !file.exists(depth_file)) {
    stop("Cannot locate depth file for summary row ", i)
  }

  depth_data <- read.delim(
    depth_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_depth_columns <- c("position", "depth")
  missing_depth_columns <- setdiff(required_depth_columns, names(depth_data))
  if (length(missing_depth_columns) > 0) {
    stop(
      "Depth file is missing column(s): ",
      paste(missing_depth_columns, collapse = ", "),
      "\nFile: ", depth_file
    )
  }
  depth_data$position <- as.numeric(depth_data$position)
  depth_data$depth <- as.numeric(depth_data$depth)
  depth_data <- depth_data[
    is.finite(depth_data$position) & is.finite(depth_data$depth),
    ,
    drop = FALSE
  ]
  depth_data <- depth_data[order(depth_data$position), , drop = FALSE]

  sequence_match <- targets$sequence_id[
    match(as.character(summary_data$bam_reference[i]), targets$bam_reference)
  ]
  if (is.na(sequence_match)) {
    stop(
      "Cannot map bam_reference to sequence_id for panel ", i, ": ",
      summary_data$bam_reference[i]
    )
  }

  panel_annotations <- annotation_data %>%
    filter(sequence_id == sequence_match)

  plot_list[[i]] <- make_depth_annotation_plot(
    depth_data = depth_data,
    species_name = summary_data$species_name[i],
    library_id = summary_data$library[i],
    reads = summary_data$reads[i],
    rpm = summary_data$rpm[i],
    coverage_pct = 100 * as.numeric(summary_data$breadth_1x[i]),
    panel_label = LETTERS[i],
    sequence_id = sequence_match,
    panel_annotations = panel_annotations
  )
}

# Render without patchwork/gridExtra. Every panel remains a regular ggplot,
# while the lower part of each panel is the gggenes protein track.
render_combined <- function(filename, device_fun, width, height, res = NULL) {
  if (is.null(res)) {
    device_fun(filename, width = width, height = height, family = "sans")
  } else {
    # For grDevices::png(), width/height default to pixels. Explicitly use
    # inches so width = 10 and height = 8.05 remain publication dimensions.
    device_fun(
      filename,
      width = width,
      height = height,
      units = "in",
      res = res
    )
  }
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 3,
    ncol = 2,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(1, 1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_len(6)) {
    if (i <= length(plot_list)) {
      print(
        plot_list[[i]],
        vp = grid::viewport(
          layout.pos.row = ceiling(i / 2),
          layout.pos.col = ifelse(i %% 2 == 1, 1, 2)
        )
      )
    }
  }
  grid::popViewport()
  grDevices::dev.off()
}

pdf_file <- file.path(outdir, "FigureS3_brain_depth_with_NR_annotations.pdf")
png_file <- file.path(outdir, "FigureS3_brain_depth_with_NR_annotations.png")

render_combined(pdf_file, grDevices::cairo_pdf, width = 10, height = 8.05)
render_combined(png_file, grDevices::png, width = 10, height = 8.05, res = 600)

message("FigureS3 depth + gggenes annotation visualization completed: ", outdir)
message("PDF: ", pdf_file)
message("PNG: ", png_file)
