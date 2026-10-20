#!/usr/bin/env Rscript

# Figure 4A: classification (inner ring) and species novelty (outer ring).
# The layout and palette follow the target Figure 4A panel supplied by the user.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggnewscale)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

input_file <- file.path(script_dir, "figure4A_kingdom_novelty_counts_input.csv")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

dat <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE,
                check.names = FALSE)

required_cols <- c("kingdom", "total_count", "known", "novel")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat$kingdom_clean <- sub("^k__", "", as.character(dat$kingdom))
dat$total_count <- as.numeric(dat$total_count)
dat$known <- as.numeric(dat$known)
dat$novel <- as.numeric(dat$novel)
if (anyNA(dat$total_count) || anyNA(dat$known) || anyNA(dat$novel)) {
  stop("Figure 4A input contains non-numeric count values.")
}

# Match the target panel orientation: Bacteria is the final horizontal
# sector on the right, Eukaryota is immediately above it, and Virus occupies
# the remaining large sector.
available_order <- c("Bacteria", "Eukaryota", "Viruses")
kingdom_order <- available_order[available_order %in% dat$kingdom_clean]
kingdom_order <- c(kingdom_order,
                   sort(setdiff(unique(dat$kingdom_clean), available_order)))

wide <- data.frame(
  kingdom = dat$kingdom_clean,
  n.known = dat$known,
  n.novel = dat$novel,
  total = dat$total_count,
  stringsAsFactors = FALSE
)
wide <- wide[match(kingdom_order, wide$kingdom), , drop = FALSE]
if (anyNA(wide$kingdom)) stop("Missing one or more expected kingdoms in Figure 4A input.")

display_name <- c(Viruses = "Virus", Eukaryota = "Eukaryota", Bacteria = "Bacteria")
wide$classification <- unname(display_name[wide$kingdom])

# Exact palette sampled from the target Figure 4A panel.
classification_colors <- c(
  Virus = "#F7B1B7",
  Eukaryota = "#C2D5B4",
  Bacteria = "#B4B4D7"
)
novelty_colors <- c(
  Known = "#E4ECEC",
  Novel = "#C5EEF7"
)

# Inner ring: classification totals.
inner <- wide
inner$ymin <- c(0, head(cumsum(inner$total), -1))
inner$ymax <- cumsum(inner$total)
inner$ymid <- (inner$ymin + inner$ymax) / 2
inner$inner_label <- paste0("n=", inner$total)
inner$classification <- factor(
  inner$classification,
  levels = c("Bacteria", "Eukaryota", "Virus")
)

# Outer ring: Known/Novel is calculated separately within each classification.
# The order must match the inner ring: Virus, Eukaryota, Bacteria.
outer <- rbind(
  data.frame(kingdom = wide$kingdom, classification = wide$classification,
             novelty = "Known", n = wide$n.known),
  data.frame(kingdom = wide$kingdom, classification = wide$classification,
             novelty = "Novel", n = wide$n.novel)
)
outer$kingdom <- factor(outer$kingdom, levels = kingdom_order)
outer$novelty <- factor(outer$novelty, levels = c("Known", "Novel"))
outer <- outer[order(outer$kingdom, outer$novelty), , drop = FALSE]
outer <- outer[outer$n > 0, , drop = FALSE]
outer$kingdom <- as.character(outer$kingdom)
outer$classification <- as.character(outer$classification)
outer$novelty <- as.character(outer$novelty)
outer$ymin <- c(0, head(cumsum(outer$n), -1))
outer$ymax <- cumsum(outer$n)
outer$ymid <- (outer$ymin + outer$ymax) / 2
outer$pct <- 100 * outer$n / wide$total[match(outer$kingdom, wide$kingdom)]
outer$outer_label <- paste0(outer$novelty, " (", sprintf("%.1f", outer$pct), "%)")
outer$novelty <- factor(outer$novelty, levels = c("Known", "Novel"))

stats_out <- data.frame(
  Classification = wide$classification,
  Known = wide$n.known,
  Novel = wide$n.novel,
  Total = wide$total,
  stringsAsFactors = FALSE
)
write.csv(
  stats_out,
  file.path(script_dir, "01.species_novelty_stats_R.csv"),
  row.names = FALSE,
  quote = FALSE
)

plot_obj <- ggplot() +
  geom_rect(
    data = inner,
    aes(xmin = 0, xmax = 0.65, ymin = ymin, ymax = ymax, fill = classification),
    color = "white",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    name = "Classification",
    values = classification_colors,
    breaks = c("Virus", "Bacteria", "Eukaryota"),
    guide = guide_legend(order = 2, nrow = 1, title.position = "top")
  ) +
  geom_text(
    data = inner,
    aes(x = 0.325, y = ymid, label = inner_label),
    size = 4,
    fontface = "bold"
  ) +
  ggnewscale::new_scale_fill() +
  geom_rect(
    data = outer,
    aes(xmin = 0.65, xmax = 1, ymin = ymin, ymax = ymax, fill = novelty),
    color = "white",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    name = "Species novelty",
    values = novelty_colors,
    breaks = c("Known", "Novel"),
    guide = guide_legend(order = 1, nrow = 1, title.position = "top")
  ) +
  geom_text(
    data = outer,
    aes(x = 1.08, y = ymid, label = outer_label),
    size = 3.2,
    hjust = 0
  ) +
  # Rotate 180 degrees to match the target panel: Bacteria at the
  # right-horizontal position, Eukaryota above it, and Virus as the main sector.
  coord_polar(theta = "y", start = 3 * pi / 2, direction = -1, clip = "off") +
  xlim(0, 1.3) +
  labs(tag = "A") +
  theme_void(base_family = "Arial") +
  theme(
    plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
    plot.tag.position = c(0.01, 0.99),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = grid::unit(0.35, "cm"),
    legend.spacing.x = grid::unit(0.35, "cm"),
    plot.margin = margin(8, 45, 8, 18)
  )

ggsave(
  file.path(script_dir, "01.species_novelty_distribution_R.pdf"),
  plot_obj,
  width = 6.5,
  height = 5.5,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  file.path(script_dir, "01.species_novelty_distribution_R.png"),
  plot_obj,
  width = 6.5,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

cat("\n=== Figure 4A summary ===\n")
print(stats_out, row.names = FALSE)
cat("Total species:", sum(wide$total), "\n")
cat("Known:", sum(wide$n.known), " Novel:", sum(wide$n.novel), "\n")
cat("Outputs saved in:", script_dir, "\n")
