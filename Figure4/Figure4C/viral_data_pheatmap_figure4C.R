# Required packages
library(openxlsx)
library(pheatmap)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
setwd(script_dir)

# Read and order the matrix
data <- read.xlsx("viral_data_pheatmap-final-input.xlsx", sheet = 1)

rownames(data) <- data[, 1]
data <- data[, -1]

desired_order <- scan("order_data.txt", what = "", sep = "\t")
valid_columns <- desired_order[desired_order %in% colnames(data)]

invalid_columns <- desired_order[!desired_order %in% colnames(data)]
if (length(invalid_columns) > 0) {
  warning("Missing matrix columns: ", paste(invalid_columns, collapse = ", "))
}
data <- data[, valid_columns]

desired_order_row <- scan("order_data_row.txt", what = "", sep = "\t")
data <- data[desired_order_row, ]

write.table(data, "data.tsv", sep = "\t", quote = FALSE)
heatmap_data <- log10(data + 1)

# Column annotations
annotation_col <- read.table("col_heatmap_annot.txt", header = TRUE, sep = "\t", row.names = 1,
                             check.names = FALSE)
annotation_col <- annotation_col[, c("Type", "Location", "Domestic_or_Wild"), drop = FALSE]
colnames(annotation_col) <- c("Sample type", "Location", "Domestic or wild")

color_col <- list(
  `Domestic or wild` = c(
    "Domestic" = "#B9A8C2",
    "Wild" = "#B0BBCB"
  ),
  Location = c(
    "Antarctica" = "#F57C6E",
    "Australia" = "#F2B56F",
    "Brazil" = "#FAE69E",
    "China" = "#84C3B7",
    "Germany" = "#88D8DB",
    "New Zealand" = "#71B7ED",
    "the United States" = "#B8AEEB"
  ),
  `Sample type` = c(
    "Cloacal swab" = "#F44236", "Mixed type libs" = "#F5CCCC", "Brain" = "#BEBADC",
    "Liver" = "#FB8072", "Lung" = "#80B1D3", "Intestine" = "#FDB462",
    "Fecal samples" = "#DC917B", "Cecal content" = "#FCCDE5", "Heart" = "#A64D79",
    "Spleen" = "#BC80BD", "Kidney" = "#CCEBC5", "Oral swab" = "#F8DE48",
    "Feathers" = "#A5CEE3", "Blood" = "#1F78B4", "Lesion" = "#B2DF8A"
  )
)

# Row annotations
annotation_row <- read.table("row_heatmap_annot.txt", header = TRUE, sep = "\t", row.names = 1,
                             check.names = FALSE)
colnames(annotation_row) <- c("Phylum", "Classification")
annotation_row$Phylum[annotation_row$Phylum == "unclassified Viruses phylum"] <- "unclassified phylum"

color_annotation_row <- list(
  Classification = c(
    "RNA_virus" = "#FBB4AE",
    "DNA_virus" = "#FFD9A8",
    "Bacteria" = "#B3CDE4",
    "Eukaryota" = "#CCEBC5"
  ),
  Phylum = c(
    "Pisuviricota" = "#CB181D", "Duplornaviricota" = "#EF3B2C", "Kitrinoviricota" = "#FB6A4A",
    "Negarnaviricota" = "#FC9272", "Artverviricota" = "#08519C", "Cossaviricota" = "#2171B5",
    "Cressdnaviricota" = "#4292C6", "Peploviricota" = "#6BAED6", "Preplasmiviricota" = "#9ECAE1",
    "unclassified phylum" = "#A64D79", "Bacillota" = "#88419D", "Bacteroidota" = "#8C6BB1",
    "Campylobacterota" = "#8C96C6", "Mycoplasmatota" = "#9EBBDA", "Pseudomonadota" = "#BFD3E6",
    "Apicomplexa" = "#004526", "Arthropoda" = "#006836", "Euglenozoa" = "#238443",
    "Microsporidia" = "#41AB5D", "Nematoda" = "#78C679", "Parabasalia" = "#ADF08E",
    "Platyhelminthes" = "#D9F3A3"
  )
)

annotation_row$Classification <- factor(annotation_row$Classification, levels = names(color_annotation_row$Classification))
annotation_row$Phylum <- factor(annotation_row$Phylum, levels = names(color_annotation_row$Phylum))

annotation_col$`Domestic or wild` <- factor(annotation_col$`Domestic or wild`, levels = names(color_col$`Domestic or wild`))
annotation_col$Location <- factor(annotation_col$Location, levels = names(color_col$Location))
annotation_col$`Sample type` <- factor(annotation_col$`Sample type`, levels = names(color_col$`Sample type`))

# Plot the heatmap
display_row_labels <- gsub("_", " ", rownames(heatmap_data), fixed = TRUE)
family_labels <- c(
  "Astroviridae", "Caliciviridae", "Hepeviridae", "Picornaviridae",
  "Sedoreoviridae", "Adenoviridae", "Parvoviridae", "Spinareoviridae",
  "Coronaviridae", "Eimeriidae", "Flaviviridae", "Plasmodiidae",
  "Circoviridae", "Paramyxoviridae", "Phenuiviridae", "Hepadnaviridae",
  "Rhabdoviridae", "Trichomonadidae", "Taeniidae", "Anelloviridae",
  "Orthoherpesviridae", "Orthomyxoviridae", "Papillomaviridae", "Kolmioviridae",
  "Pneumoviridae", "Campylobacteraceae", "Clostridiaceae", "Mycoplasmataceae",
  "Pasteurellaceae", "Staphylococcaceae", "Streptococcaceae", "Vibrionaceae",
  "Weeksellaceae", "Acaridae", "Cryptosporidiidae", "Hymenolepididae",
  "Pleistophoridae", "Sarcoptidae", "Strongyloididae", "Tritrichomonadidae",
  "Trombiculidae", "Trypanosomatidae"
)

df <- pheatmap(
  heatmap_data, 
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = list(
    `Sample type` = color_col$`Sample type`,
    Location = color_col$Location,
    `Domestic or wild` = color_col$`Domestic or wild`,
    Phylum = color_annotation_row$Phylum,
    Classification = color_annotation_row$Classification
  ), 
  color = colorRampPalette(c("#f7f6f5", "#DD1332"))(200),
  show_colnames = FALSE,
  labels_row = display_row_labels,
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  gaps_row = c(9, 17, 24, 40, 44, 52, 57, 58, 59, 61, 64, 65, 66, 67, 69, 82, 85, 86, 87, 88, 92),
  gaps_col = c(2, 41, 42, 245, 249, 354),
  fontsize_row = 11,
  legend = TRUE,
  legend_position = "bottom"
)

# Format taxon labels.
row_name_idx <- which(df$layout$name == "row_names")
if (length(row_name_idx) == 1) {
  row_name_grob <- df$grobs[[row_name_idx]]
  row_name_grob$gp <- grid::gpar(
    fontsize = 11,
    fontface = ifelse(display_row_labels %in% family_labels, "bold.italic", "italic"),
    col = "black"
  )
  df$grobs[[row_name_idx]] <- row_name_grob
}

row_annotation_name_idx <- grep("row_annotation_names|annotation_names_row", df$layout$name)
for (i in row_annotation_name_idx) {
  if (inherits(df$grobs[[i]], "text")) {
    df$grobs[[i]]$rot <- 0
    df$grobs[[i]]$hjust <- 0.5
    df$grobs[[i]]$vjust <- 0.5
  }
}

pdf_file <- file.path(getwd(), "heatmap_output_figure4C.pdf")
pdf(pdf_file, width = 18, height = 21)
print(df)
dev.off()

message("Saved: ", pdf_file)
