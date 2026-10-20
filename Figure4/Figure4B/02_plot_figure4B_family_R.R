#!/usr/bin/env Rscript

# Figure 4B: number of unique species per family.
# Layout follows the supplied target panel: horizontal bars, log10 x-axis,
# italic family labels, counts at bar ends, and a classification color strip.

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

input_file <- file.path(script_dir, "figure4B_family_counts_input.csv")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

family_counts <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("rank", "family", "count", "classification")
missing_cols <- setdiff(required_cols, names(family_counts))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

family_counts$rank <- as.integer(family_counts$rank)
family_counts$species_count <- as.numeric(family_counts$count)
family_counts$family <- sub("^f__", "", as.character(family_counts$family))
family_counts$classification <- as.character(family_counts$classification)
family_counts <- family_counts[order(family_counts$rank), , drop = FALSE]

if (nrow(family_counts) == 0) stop("No family/species records available after filtering.")

# Highest-count families appear at the top of the chart.
family_order <- family_counts$family
family_counts$family <- factor(family_counts$family, levels = rev(family_order))
family_counts$y_pos <- nrow(family_counts):1

# Target Figure 4B palette.
bar_color <- "#68689B"
classification_colors <- c(
  Virus = "#F7B5B9",
  Eukaryota = "#C8DCC4",
  Bacteria = "#ADB1D1"
)

# The strip is positioned to the right of the bars on the log scale.
# Use explicit rectangles rather than geom_tile so that the annotation is
# retained reliably when the y axis is discrete and the figure is exported.
strip_xmin <- 360
strip_xmax <- 470
strip_label_x <- 610
strip_y <- nrow(family_counts):1
classification_strip <- data.frame(
  classification = family_counts$classification,
  xmin = strip_xmin,
  xmax = strip_xmax,
  ymin = strip_y - 0.39,
  ymax = strip_y + 0.39,
  stringsAsFactors = FALSE
)

write.csv(
  family_counts[, c("rank", "family", "classification", "species_count")],
  file.path(script_dir, "figure4B_family_species_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

plot_obj <- ggplot(family_counts, aes(y = y_pos)) +
  geom_col(
    aes(x = species_count),
    width = 0.78,
    fill = bar_color,
    orientation = "y"
  ) +
  geom_text(
    aes(x = species_count * 1.04, label = species_count),
    hjust = 0,
    size = 3.1
  ) +
  geom_rect(
    data = classification_strip,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = classification),
    inherit.aes = FALSE,
    color = NA
  ) +
  scale_fill_manual(values = classification_colors, guide = "none") +
  annotate(
    "text",
    x = strip_label_x,
    y = mean(seq_len(nrow(family_counts))),
    label = "Classification",
    angle = 90,
    fontface = "bold",
    size = 4.2
  ) +
  scale_x_log10(
    breaks = c(1, 10, 100),
    labels = c("1", "10", "100"),
    limits = c(0.8, 700),
    position = "top",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = family_counts$y_pos,
    labels = as.character(family_counts$family),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  labs(tag = "B", x = "Species count", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.tag = element_text(size = 17, face = "bold", hjust = 0, vjust = 1),
    plot.tag.position = c(0.01, 0.99),
    axis.text.y = element_text(size = 8.5, face = "italic", color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    plot.margin = margin(8, 55, 12, 12)
  )

ggsave(
  file.path(script_dir, "Figure4B_family_species_count_R.pdf"),
  plot_obj,
  width = 5.5,
  height = max(8, 0.28 * nrow(family_counts) + 2.5),
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  file.path(script_dir, "Figure4B_family_species_count_R.png"),
  plot_obj,
  width = 5.5,
  height = max(8, 0.28 * nrow(family_counts) + 2.5),
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

cat("\n=== Figure 4B summary ===\n")
print(family_counts[, c("rank", "family", "species_count", "classification")], row.names = FALSE)
cat("Families plotted:", nrow(family_counts), "\n")
cat("Outputs saved in:", script_dir, "\n")
