#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(readxl)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
setwd(script_dir)

clean <- function(x) {
  y <- as.character(x)
  y[is.na(x)] <- ""
  trimws(y)
}

rpm <- as.data.frame(readxl::read_excel(
  "Figure5A_final_Table5_RPM.xlsx", sheet = "main", skip = 1,
  .name_repair = "unique"
), stringsAsFactors = FALSE)
annotation <- as.data.frame(readxl::read_excel(
  "Figure5A_final_Table4_annotation.xlsx", sheet = "Sheet2", skip = 1,
  .name_repair = "unique"
), stringsAsFactors = FALSE)
metadata <- as.data.frame(readxl::read_excel(
  "Figure5A_final_Table3_metadata.xlsx", sheet = "Sheet4", skip = 1,
  .name_repair = "unique"
), stringsAsFactors = FALSE)

names(rpm) <- trimws(names(rpm))
names(annotation) <- trimws(names(annotation))
names(metadata) <- trimws(names(metadata))

required_rpm <- c("lib_id", "kindom", "family", "species_name")
required_annotation <- c("lib_contig", "phylum")
required_metadata <- c("SRA/GSA RUN NUMBER", "Sample Type")
if (!all(required_rpm %in% names(rpm))) stop("Table5 is missing required columns.")
if (!all(required_annotation %in% names(annotation))) stop("Table4 is missing required columns.")
if (!all(required_metadata %in% names(metadata))) stop("Table3 is missing required columns.")

rpm$lib_id <- clean(rpm$lib_id)
annotation$lib_contig <- clean(annotation$lib_contig)
metadata$lib_id <- clean(metadata[["SRA/GSA RUN NUMBER"]])
metadata$original_sample_type <- clean(metadata[["Sample Type"]])

if (anyDuplicated(rpm$lib_id)) stop("Table5 contains duplicated reference IDs.")
if (anyDuplicated(annotation$lib_contig)) stop("Table4 contains duplicated reference IDs.")
if (anyDuplicated(metadata$lib_id)) stop("Table3 contains duplicated library IDs.")

library_columns <- setdiff(names(rpm), c("No.", required_rpm))
if (!setequal(library_columns, metadata$lib_id)) {
  stop("Table5 library IDs and Table3 library IDs do not match.")
}
if (!setequal(rpm$lib_id, annotation$lib_contig)) {
  stop("Table4 and Table5 reference IDs do not match.")
}

sample_map <- c(
  "Cloacal swab" = "Cloacal swab",
  "Oral swab" = "Oral swab",
  "Fecal samples" = "Fecal samples",
  "Stool samples" = "Fecal samples",
  "feces" = "Fecal samples",
  "Cecal content" = "Cecal content",
  "Cecal contents" = "Cecal content",
  "Cecal samples" = "Cecal content",
  "Intestine" = "Intestine",
  "Intestinal samples" = "Intestine",
  "gut" = "Intestine",
  "Liver" = "Liver",
  "Lung" = "Lung",
  "Brain" = "Brain",
  "Spleen" = "Spleen",
  "Kidney" = "Kidney",
  "Heart" = "Heart",
  "Blood" = "Blood",
  "blood" = "Blood",
  "Lesion" = "Lesion",
  "Dorsal feathers" = "Feathers",
  "Ventral feathers" = "Feathers"
)
metadata$sample_type <- unname(sample_map[metadata$original_sample_type])

numeric <- as.data.frame(lapply(rpm[library_columns], function(x) {
  z <- suppressWarnings(as.numeric(x))
  if (anyNA(z) || any(z < 0)) stop("RPM matrix contains invalid values.")
  z[z < 1] <- 0
  z
}), check.names = FALSE, stringsAsFactors = FALSE)

sample_order <- c(
  "Cloacal swab", "Oral swab", "Fecal samples", "Cecal content",
  "Liver", "Intestine", "Brain", "Lung", "Spleen", "Kidney", "Heart",
  "Blood", "Lesion"
)

sample_matrix <- do.call(cbind, lapply(sample_order, function(sample_type) {
  ids <- metadata$lib_id[!is.na(metadata$sample_type) & metadata$sample_type == sample_type]
  if (!length(ids)) return(rep(0, nrow(rpm)))
  rowSums(numeric[, ids, drop = FALSE])
}))
colnames(sample_matrix) <- sample_order

annotation_map <- annotation$phylum[match(rpm$lib_id, annotation$lib_contig)]
annotation_map <- clean(annotation_map)
annotation_map[annotation_map == ""] <- "unclassified viruses"

matrix_data <- data.frame(
  lib_id = rpm$lib_id,
  kindom = clean(rpm$kindom),
  phylum = annotation_map,
  family = clean(rpm$family),
  species_name = clean(rpm$species_name),
  sample_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

abundance <- as.matrix(matrix_data[, sample_order, drop = FALSE])
storage.mode(abundance) <- "numeric"
phylum_sum <- rowsum(abundance, matrix_data$phylum, reorder = FALSE)
phylum_total <- colSums(phylum_sum)
phylum_relative <- sweep(phylum_sum, 2, phylum_total, "/")
phylum_relative[is.na(phylum_relative)] <- 0

infection_counts <- colSums(abundance > 0)
counts <- as.numeric(infection_counts)
names(counts) <- names(infection_counts)

phylum_colors <- c(
  Apicomplexa = "#FDDFE8", Arthropoda = "#E67D90", Artverviricota = "#D92F3D",
  Bacillota = "#F39FC9", Bacteroidota = "#DF2C79", Campylobacterota = "#C91D63",
  Cossaviricota = "#9E1E5B", Cressdnaviricota = "#D5A9DA", Duplornaviricota = "#B87DC5",
  Euglenozoa = "#9C57B2", Kitrinoviricota = "#8941A5", Microsporidia = "#5D2B95",
  Mycoplasmatota = "#472C85", Negarnaviricota = "#DCEFFA", Nematoda = "#B1D8F3",
  Parabasalia = "#4C8CD1", Peploviricota = "#2D4CA1", Pisuviricota = "#A9E2E2",
  Platyhelminthes = "#52C5CB", Preplasmiviricota = "#2C5E66",
  Pseudomonadota = "#A8C137", `unclassified viruses` = "#D2D0E6"
)

phylum_labels <- sub("^p__", "", rownames(phylum_relative))
phylum_labels[phylum_labels == "unclassified Viruses phylum"] <- "unclassified viruses"
missing_phyla <- setdiff(phylum_labels, names(phylum_colors))
if (length(missing_phyla) > 0) phylum_colors[missing_phyla] <- "#BDBDBD"

relative_long <- do.call(rbind, lapply(seq_along(sample_order), function(j) {
  data.frame(
    Sample_Type = sample_order[j],
    phylum = phylum_labels,
    relative = phylum_relative[, j],
    stringsAsFactors = FALSE
  )
}))
relative_long$Sample_Type <- factor(relative_long$Sample_Type, levels = sample_order)
relative_long$x <- match(as.character(relative_long$Sample_Type), sample_order)
relative_long$phylum <- factor(relative_long$phylum, levels = names(phylum_colors))

group_data <- data.frame(
  group = c("Swabs", "Luminal contents", "Organs", "Fluids", "Lesion"),
  xmin = c(0.55, 2.55, 4.55, 11.55, 12.55),
  xmax = c(2.45, 4.45, 11.45, 12.45, 13.45),
  stringsAsFactors = FALSE
)

sample_labels <- data.frame(
  x = seq_along(sample_order) - 0.34,
  label = sample_order,
  stringsAsFactors = FALSE
)

panel_xmin <- 0.5
panel_xmax <- length(sample_order) + 0.5

bar_border <- data.frame(
  x = c(panel_xmin, panel_xmax, panel_xmin, panel_xmin),
  xend = c(panel_xmax, panel_xmax, panel_xmin, panel_xmax),
  y = c(0, 0, 0, 1),
  yend = c(0, 1, 1, 1)
)

top_plot <- ggplot(relative_long, aes(x = x, y = relative, fill = phylum)) +
  geom_col(width = 0.78, colour = NA) +
  geom_segment(data = bar_border,
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.35, colour = "black") +
  geom_text(data = sample_labels,
            aes(x = x, y = 1.015, label = label),
            inherit.aes = FALSE, angle = 45, hjust = 0, vjust = 0,
            size = 3.2) +
  geom_segment(data = group_data,
               aes(x = xmin, xend = xmax, y = 1.30, yend = 1.30),
               inherit.aes = FALSE, linewidth = 0.35) +
  geom_text(data = group_data,
            aes(x = (xmin + xmax) / 2, y = 1.345, label = group),
            inherit.aes = FALSE, fontface = "bold", size = 3.8) +
  annotate("text", x = 0.58, y = 1.40, label = "A", fontface = "bold", size = 6) +
  scale_fill_manual(values = phylum_colors, breaks = names(phylum_colors), name = "Phylum") +
  scale_x_continuous(breaks = NULL, limits = c(panel_xmin, panel_xmax)) +
  scale_y_continuous(limits = c(0, 1.45), breaks = seq(0, 1, 0.2),
                     labels = function(x) format(x, trim = TRUE), expand = c(0, 0)) +
  labs(x = NULL, y = "Relative RPM percentage (%)") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"), legend.text = element_text(size = 8),
        legend.key.height = grid::unit(0.34, "cm"),
        legend.key.width = grid::unit(0.28, "cm"),
        plot.margin = margin(5, 0, 0, 0)) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

count_data <- data.frame(x = seq_along(sample_order), count = counts[sample_order])
count_data$label_y <- -pmax(count_data$count / 2, 7)
count_max <- ceiling(max(count_data$count, na.rm = TRUE) / 50) * 50

bottom_plot <- ggplot(count_data, aes(x = x, y = -count)) +
  geom_col(width = 0.78, fill = "#AFAFAF") +
  geom_text(aes(y = label_y, label = count), size = 3.4, colour = "black") +
  scale_x_continuous(breaks = seq_along(sample_order), labels = NULL,
                     limits = c(panel_xmin, panel_xmax)) +
  scale_y_continuous(limits = c(-count_max * 1.08, 0),
                     breaks = seq(0, -count_max, by = -50),
                     labels = function(x) abs(x), expand = c(0, 0)) +
  labs(x = NULL, y = "Species count") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_blank(), plot.margin = margin(0, 0, 5, 0))

top_core <- top_plot + theme(legend.position = "none")
legend_grob <- gtable::gtable_filter(ggplotGrob(top_plot), "guide-box")
top_grob <- ggplotGrob(top_core)
bottom_grob <- ggplotGrob(bottom_plot)
shared_widths <- grid::unit.pmax(top_grob$widths, bottom_grob$widths)
top_grob$widths <- shared_widths
bottom_grob$widths <- shared_widths
main_panel <- arrangeGrob(top_grob, grid::nullGrob(), bottom_grob,
                          ncol = 1, heights = c(2.1, 0.12, 1.55))
figure <- arrangeGrob(main_panel, legend_grob, ncol = 2, widths = c(1, 0.22))

pdf("Figure5A_tissue_landscape.pdf", width = 12, height = 7.2, useDingbats = FALSE)
grid::grid.draw(figure)
dev.off()
ggsave("Figure5A_tissue_landscape.png", figure, width = 12, height = 7.2,
       units = "in", dpi = 300)
