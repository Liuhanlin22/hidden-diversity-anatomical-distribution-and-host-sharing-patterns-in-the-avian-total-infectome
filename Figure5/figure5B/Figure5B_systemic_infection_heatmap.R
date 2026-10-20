suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
setwd(script_dir)

unlink(c(
  ".RData",
  ".Rhistory",
  "classification_corrections.csv",
  "classification_counts.csv",
  "virus_table_systemic_tissues_by_species.csv",
  "virus_table_systemic_tissues_final.csv"
), force = TRUE)

metadata_file <- "Figure5B_final_Table3_metadata.xlsx"
rpm_file <- "Figure5B_final_Table5_RPM.xlsx"
output_pdf <- "Figure5B_systemic_infection_heatmap.pdf"
output_png <- "Figure5B_systemic_infection_heatmap.png"

if (!file.exists(metadata_file)) stop("Input file not found: ", metadata_file)
if (!file.exists(rpm_file)) stop("Input file not found: ", rpm_file)

metadata <- read_excel(metadata_file, sheet = "Sheet4", skip = 1)
rpm <- read_excel(rpm_file, sheet = "main", skip = 1)

required_metadata <- c("SRA/GSA RUN NUMBER", "Sample Type")
required_rpm <- c("No.", "lib_id", "kindom", "family", "species_name")
if (!all(required_metadata %in% names(metadata))) {
  stop("Metadata columns are incomplete.")
}
if (!all(required_rpm %in% names(rpm))) {
  stop("RPM matrix columns are incomplete.")
}

library_ids <- names(rpm)[6:ncol(rpm)]
metadata_ids <- trimws(as.character(metadata[["SRA/GSA RUN NUMBER"]]))
metadata_types <- trimws(as.character(metadata[["Sample Type"]]))
metadata_map <- setNames(metadata_types, metadata_ids)
tissue_by_library <- unname(metadata_map[library_ids])
tissues <- c("Brain", "Heart", "Lung", "Kidney", "Liver", "Spleen", "Intestine", "Blood")
if (any(is.na(tissue_by_library))) {
  missing_metadata <- library_ids[is.na(tissue_by_library)]
  stop("Missing metadata for libraries: ", paste(missing_metadata, collapse = ", "))
}

normalize_tissue <- function(x) {
  x <- tolower(trimws(x))
  x[x %in% c("brain")] <- "Brain"
  x[x %in% c("heart")] <- "Heart"
  x[x %in% c("lung")] <- "Lung"
  x[x %in% c("kidney")] <- "Kidney"
  x[x %in% c("liver")] <- "Liver"
  x[x %in% c("spleen")] <- "Spleen"
  x[x %in% c("intestine", "intestinal samples", "gut")] <- "Intestine"
  x[x %in% c("blood")] <- "Blood"
  x[!x %in% tissues] <- NA_character_
  x
}
tissue_by_library <- normalize_tissue(tissue_by_library)

rpm_matrix <- as.matrix(as.data.frame(lapply(rpm[library_ids], function(x) {
  as.numeric(as.character(x))
})))
storage.mode(rpm_matrix) <- "numeric"
rpm_matrix[is.na(rpm_matrix)] <- 0
rpm_matrix[rpm_matrix < 1] <- 0

tissue_values <- sapply(tissues, function(tissue) {
  selected <- which(tissue_by_library == tissue)
  if (length(selected) == 0) {
    rep(0, nrow(rpm_matrix))
  } else {
    rowSums(rpm_matrix[, selected, drop = FALSE], na.rm = TRUE)
  }
})
tissue_values <- as.data.frame(tissue_values, check.names = FALSE)

data <- data.frame(
  classification = as.character(rpm$kindom),
  family = as.character(rpm$family),
  species_name = as.character(rpm$species_name),
  tissue_values,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

selected_rows <- apply(data[tissues], 1, function(x) any(x > 0))
data <- data[selected_rows, , drop = FALSE]
if (nrow(data) == 0) stop("No species have nonzero RPM in the selected tissues.")
if (anyDuplicated(data$species_name)) stop("Duplicate species_name values detected.")

display_name_map <- c(
  "Alphainfluenzavirus influenzae (BJ94-Like H9)" = "BJ94-Like H9 influenza A virus",
  "Alphainfluenzavirus influenzae (ST2853-Like H6)" = "ST2853-Like H6 influenza A virus",
  "Alphainfluenzavirus influenzae (2.3.4.4h H5)" = "2.3.4.4h H5 influenza A virus",
  "Avian associated chaphamaparvovirus sp. L1" = "Avian associated chaphamaparvo- sp. L1",
  "Avian associated chaphamaparvovirus sp. L2" = "Avian associated chaphamaparvo- sp. L2",
  "Avian associated chaphamaparvovirus sp. L3" = "Avian associated chaphamaparvo- sp. L3",
  "Avian associated chaphamaparvovirus sp. L4" = "Avian associated chaphamaparvo- sp. L4",
  "Avian associated chaphamaparvovirus sp. L5" = "Avian associated chaphamaparvo- sp. L5",
  "Avian associated chaphamaparvovirus sp. L6" = "Avian associated chaphamaparvo- sp. L6"
)
mapped_names <- unname(display_name_map[data$species_name])
data$species_name[!is.na(mapped_names)] <- mapped_names[!is.na(mapped_names)]

data$family_label <- sub("^f__", "", data$family)
family_class_map <- c(
  Astroviridae = "RNA_virus",
  Caliciviridae = "RNA_virus",
  Flaviviridae = "RNA_virus",
  Hepeviridae = "RNA_virus",
  Orthomyxoviridae = "RNA_virus",
  Paramyxoviridae = "RNA_virus",
  Picornaviridae = "RNA_virus",
  Rhabdoviridae = "RNA_virus",
  Adenoviridae = "DNA_virus",
  Anelloviridae = "DNA_virus",
  Circoviridae = "DNA_virus",
  Hepadnaviridae = "DNA_virus",
  Parvoviridae = "DNA_virus",
  Eimeriidae = "k__Eukaryota",
  Taeniidae = "k__Eukaryota",
  Plasmodiidae = "k__Eukaryota",
  Tritrichomonadidae = "k__Eukaryota",
  Trypanosomatidae = "k__Eukaryota"
)
data$classification <- ifelse(
  data$classification %in% c("k__Eukaryota", "k__Bacteria"),
  data$classification,
  unname(family_class_map[data$family_label])
)
if (anyNA(data$classification)) {
  stop("Unmapped family classification: ",
       paste(unique(data$family_label[is.na(data$classification)]), collapse = ", "))
}
class_order <- c("RNA_virus", "DNA_virus", "k__Bacteria", "k__Eukaryota")
data$classification <- factor(data$classification, levels = class_order)
data <- data[order(data$classification, data$family_label, data$species_name), , drop = FALSE]
data$classification <- as.character(data$classification)

class_labels <- c(
  RNA_virus = "RNA virus",
  DNA_virus = "DNA virus",
  k__Eukaryota = "Eukaryota"
)
class_colors <- c(
  "RNA virus" = "#F6AAA5",
  "DNA virus" = "#FCD39B",
  Eukaryota = "#B9DEAE"
)
family_colors <- c(
  Astroviridae = "#9C0033",
  Caliciviridae = "#D50032",
  Flaviviridae = "#ED5564",
  Hepeviridae = "#FA4D3A",
  Orthomyxoviridae = "#FF8436",
  Paramyxoviridae = "#FFB54A",
  Picornaviridae = "#FFD36A",
  Rhabdoviridae = "#FBE69A",
  Circoviridae = "#3F007D",
  Hepadnaviridae = "#6A51A3",
  Parvoviridae = "#807DBA",
  Adenoviridae = "#9E9AC8",
  Anelloviridae = "#DADAEB",
  Eimeriidae = "#006D2C",
  Taeniidae = "#238443",
  Plasmodiidae = "#41AB5D",
  Tritrichomonadidae = "#78C679",
  Trypanosomatidae = "#C2E699"
)
novel_pathogens <- c(
  "Avian associated calicivirus sp. 19",
  "Avian hepacivirus sp. 1",
  "Avian hepacivirus sp. 2",
  "Avian pegivirus sp. 1",
  "Avian associated metaavulavirus sp. L3",
  "Avian associated hepato-like virus sp. L211",
  "Avian associated rhabdovirus sp. L1",
  "Avian associated chaphamaparvo- sp. L2",
  "Avian associated chaphamaparvo- sp. L3",
  "Avian associated chaphamaparvo- sp. L4",
  "Avian associated chaphamaparvo- sp. L5",
  "Avian associated chaphamaparvo- sp. L6",
  "Aviadenovirus sp. L7",
  "Aviadenovirus sp. L8",
  "Hydatigera sp. Avi1",
  "Plasmodium sp. Avi1",
  "Plasmodium sp. Avi2",
  "Plasmodium sp. Avi3",
  "Plasmodium sp. Avi4",
  "Plasmodium sp. Avi5",
  "Tritrichomonas sp . Avi1"
)

data$classification_label <- unname(class_labels[data$classification])
values <- as.matrix(data[, tissues, drop = FALSE])
storage.mode(values) <- "numeric"
values <- log10(values + 1)
nr <- nrow(data)
nc <- length(tissues)
class_break <- c(data$classification[-nr] != data$classification[-1], FALSE)
gap_after <- ifelse(class_break, 0.45, 0)
row_y <- nr - seq_len(nr) + 1 - c(0, cumsum(gap_after[-nr]))
top_y <- max(row_y)
bottom_y <- min(row_y)
heat_colors <- colorRampPalette(c(
  "#F1F0EE", "#F6C4CA", "#E88A98", "#DB5B70", "#D62B46", "#C3002F"
))(101)
value_to_color <- function(x) heat_colors[floor(pmax(0, pmin(x, 5)) / 5 * 100) + 1]

heatmap_data <- expand.grid(row = seq_len(nr), col = seq_len(nc))
heatmap_data$x <- heatmap_data$col + 2.65
heatmap_data$y <- row_y[heatmap_data$row]
heatmap_data$fill <- value_to_color(values[cbind(heatmap_data$row, heatmap_data$col)])
classification_data <- data.frame(
  x = 1, y = row_y, fill = unname(class_colors[data$classification_label])
)
family_data <- data.frame(
  x = 1.82, y = row_y, fill = unname(family_colors[data$family_label])
)
row_labels <- data.frame(
  x = nc + 3.35,
  y = row_y,
  label = data$species_name,
  novel = data$species_name %in% novel_pathogens,
  stringsAsFactors = FALSE
)

make_star <- function(x, y, radius, group_id) {
  angles <- pi / 2 + seq(0, 2 * pi, length.out = 11)[-11]
  radii <- rep(c(radius, radius * 0.42), 5)
  x_radius <- radius * 0.30
  data.frame(
    x = x + x_radius * (radii / radius) * cos(angles),
    y = y + radii * sin(angles),
    group = group_id
  )
}
novel_rows <- which(row_labels$novel)
star_data <- if (length(novel_rows) > 0) {
  do.call(rbind, lapply(seq_along(novel_rows), function(i) {
    row_id <- novel_rows[i]
    star_x <- row_labels$x[row_id] + min(7.0, 0.14 * nchar(row_labels$label[row_id]) + 0.20)
    make_star(star_x, row_labels$y[row_id], 0.42, i)
  }))
} else {
  data.frame(x = numeric(0), y = numeric(0), group = numeric(0))
}

classification_legend <- data.frame(
  label = names(class_colors),
  y = c(93, 90.6, 88.2),
  fill = unname(class_colors),
  stringsAsFactors = FALSE
)
family_legend <- data.frame(
  label = names(family_colors),
  y = seq(81.5, 53.6, length.out = length(family_colors)),
  fill = unname(family_colors),
  stringsAsFactors = FALSE
)
log_legend <- data.frame(
  x = 0.34,
  y = seq(20.5, 34.5, length.out = 101),
  fill = heat_colors
)
log_ticks <- data.frame(
  x = 0.64,
  y = seq(20.5, 34.5, length.out = 6),
  label = 0:5
)
tissue_labels <- data.frame(x = seq_len(nc) + 2.65, y = top_y + 1.15, label = tissues)
organ_columns <- seq_len(nc - 1) + 2.65

legend_star <- make_star(0.34, 44.5, 0.85, 999)
legend_panel <- ggplot() +
  geom_text(aes(x = 0.05, y = 98, label = "Classification"),
            hjust = 0, fontface = "bold", size = 3.8) +
  geom_tile(data = classification_legend, aes(x = 0.34, y = y, fill = fill),
            width = 0.34, height = 1.85, show.legend = FALSE) +
  geom_text(data = classification_legend, aes(x = 0.68, y = y, label = label),
            hjust = 0, size = 3.4, colour = "#1A1A1A") +
  geom_text(aes(x = 0.05, y = 83.5, label = "Family"),
            hjust = 0, fontface = "bold", size = 3.8) +
  geom_tile(data = family_legend, aes(x = 0.34, y = y, fill = fill),
            width = 0.34, height = 1.50, show.legend = FALSE) +
  geom_text(data = family_legend, aes(x = 0.68, y = y, label = label),
            hjust = 0, family = "Arial", fontface = "italic", size = 3.0,
            colour = "#1A1A1A") +
  geom_text(aes(x = 0.05, y = 48.5, label = "Novel pathogens"),
            hjust = 0, fontface = "bold", size = 3.4) +
  geom_polygon(data = legend_star, aes(x = x, y = y, group = group),
               fill = "#EF3B2C", colour = NA) +
  geom_text(aes(x = 0.05, y = 37.5, label = "log10(RPM+1)"),
            hjust = 0, fontface = "bold", size = 3.2) +
  geom_tile(data = log_legend, aes(x = x, y = y, fill = fill),
            width = 0.34, height = 0.14, show.legend = FALSE) +
  geom_text(data = log_ticks, aes(x = x, y = y, label = label),
            hjust = 0, size = 3.0, colour = "#1A1A1A") +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 5.8), ylim = c(0, 100), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(7, 0, 7, 7))

main_plot <- ggplot() +
  geom_tile(data = heatmap_data, aes(x = x, y = y, fill = fill),
            width = 0.86, height = 1.00, show.legend = FALSE) +
  geom_tile(data = classification_data, aes(x = x, y = y, fill = fill),
            width = 0.62, height = 1.00, show.legend = FALSE) +
  geom_tile(data = family_data, aes(x = x, y = y, fill = fill),
            width = 0.62, height = 1.00, show.legend = FALSE) +
  geom_text(data = row_labels, aes(x = x, y = y, label = label, colour = novel),
            hjust = 0, family = "Arial", fontface = "italic", size = 3.2,
            show.legend = FALSE) +
  geom_polygon(data = star_data, aes(x = x, y = y, group = group),
               fill = "#EF3B2C", colour = NA) +
  geom_text(data = tissue_labels, aes(x = x, y = y, label = label),
            angle = 45, hjust = 0, vjust = 0.5, family = "Arial",
            fontface = "bold.italic", size = 3.4) +
  geom_text(aes(x = mean(organ_columns), y = top_y + 4.55, label = "Organs"),
            hjust = 0.5, fontface = "bold", size = 3.8) +
  geom_segment(aes(x = min(organ_columns) - 0.5, xend = max(organ_columns) + 0.5,
                   y = top_y + 3.55, yend = top_y + 3.55),
               linewidth = 0.45, colour = "#1A1A1A") +
  scale_fill_identity() +
  scale_colour_manual(values = c("FALSE" = "#1A1A1A", "TRUE" = "#EF3B2C")) +
  coord_cartesian(xlim = c(0.4, nc + 11.1), ylim = c(bottom_y - 1.0, top_y + 6.8), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(12, 12, 12, 0))

figure <- arrangeGrob(legend_panel, main_plot, ncol = 2, widths = c(2.35, 11.65))

if (capabilities("cairo")) {
  grDevices::cairo_pdf(output_pdf, width = 14, height = 10)
} else {
  grDevices::pdf(output_pdf, width = 14, height = 10, useDingbats = FALSE)
}
grid::grid.draw(figure)
dev.off()
ggsave(output_png, plot = figure, width = 14, height = 10, units = "in", dpi = 300, bg = "white")
