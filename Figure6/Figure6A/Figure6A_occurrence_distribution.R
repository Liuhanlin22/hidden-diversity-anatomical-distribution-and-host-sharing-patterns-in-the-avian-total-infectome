#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required packages before running: ",
       paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(ggplot2))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) return(dirname(normalizePath(active_path)))
  }
  getwd()
}

base_dir <- get_script_dir()
input_file <- file.path(base_dir, "species_occurrence_distribution.csv")
if (!file.exists(input_file)) stop("Missing input file: ", input_file)

occurrence <- read.csv(input_file, check.names = FALSE,
                       stringsAsFactors = FALSE)
required_columns <- c("Occurrence_Category", "Count_of_Species")
missing_columns <- setdiff(required_columns, names(occurrence))
if (length(missing_columns) > 0) {
  stop("Input is missing: ", paste(missing_columns, collapse = ", "))
}

category_levels <- c(
  "1 country",
  "2 countries",
  "3 countries",
  ">= 4 countries"
)
category_labels <- c(
  "Within country",
  "2 countries",
  "3 countries",
  "\u2265 4 countries"
)
category_colors <- c(
  "Within country" = "#E4E6F2",
  "2 countries" = "#B3B5D2",
  "3 countries" = "#8F90BA",
  "\u2265 4 countries" = "#7072AD"
)

occurrence$Occurrence_Category <- trimws(as.character(occurrence$Occurrence_Category))
occurrence$Count_of_Species <- suppressWarnings(
  as.numeric(occurrence$Count_of_Species)
)
if (anyNA(occurrence$Count_of_Species) ||
    any(occurrence$Count_of_Species < 0) ||
    any(occurrence$Count_of_Species %% 1 != 0)) {
  stop("Count_of_Species must contain non-negative integer values.")
}
if (anyDuplicated(occurrence$Occurrence_Category)) {
  stop("Occurrence_Category must contain one row per category.")
}
missing_categories <- setdiff(category_levels, occurrence$Occurrence_Category)
if (length(missing_categories) > 0) {
  stop("Input is missing categories: ", paste(missing_categories, collapse = ", "))
}

occurrence <- occurrence[match(category_levels, occurrence$Occurrence_Category), , drop = FALSE]
occurrence$Category <- factor(category_labels, levels = category_labels)
occurrence$Fraction <- occurrence$Count_of_Species / sum(occurrence$Count_of_Species)
occurrence$Label_y <- cumsum(occurrence$Count_of_Species) -
  occurrence$Count_of_Species / 2
occurrence$Label <- as.character(occurrence$Count_of_Species)

write.csv(
  occurrence[, c("Category", "Count_of_Species", "Fraction")],
  file.path(base_dir, "species_occurrence_distribution_final.csv"),
  row.names = FALSE
)

pie_order <- c(
  "\u2265 4 countries",
  "3 countries",
  "2 countries",
  "Within country"
)
pie_data <- occurrence[match(pie_order, occurrence$Category), , drop = FALSE]
pie_data$angle_start <- c(0, cumsum(head(pie_data$Fraction, -1))) * 2 * pi
pie_data$angle_end <- cumsum(pie_data$Fraction) * 2 * pi

make_sector <- function(i, n = 121) {
  angles <- seq(pie_data$angle_start[i], pie_data$angle_end[i], length.out = n)
  data.frame(
    Category = pie_data$Category[i],
    x = c(0, cos(angles)),
    y = c(0, sin(angles))
  )
}

pie_polygons <- do.call(
  rbind,
  lapply(seq_len(nrow(pie_data)), make_sector)
)
pie_polygons$Category <- factor(pie_polygons$Category,
                                levels = category_labels)

label_angles <- c(
  "Within country" = pi,
  "2 countries" = (
    pie_data$angle_start[pie_data$Category == "2 countries"] +
      pie_data$angle_end[pie_data$Category == "2 countries"]
  ) / 2,
  "3 countries" = (
    pie_data$angle_start[pie_data$Category == "3 countries"] +
      pie_data$angle_end[pie_data$Category == "3 countries"]
  ) / 2,
  "\u2265 4 countries" = (
    pie_data$angle_start[pie_data$Category == "\u2265 4 countries"] +
      pie_data$angle_end[pie_data$Category == "\u2265 4 countries"]
  ) / 2
)
label_radius <- c(
  "Within country" = 0.40,
  "2 countries" = 0.76,
  "3 countries" = 0.84,
  "\u2265 4 countries" = 0.91
)
label_data <- occurrence[, c("Category", "Label"), drop = FALSE]
label_data$angle <- unname(label_angles[as.character(label_data$Category)])
label_data$radius <- unname(label_radius[as.character(label_data$Category)])
label_data$x <- label_data$radius * cos(label_data$angle)
label_data$y <- label_data$radius * sin(label_data$angle)

occurrence_plot <- ggplot() +
  geom_polygon(data = pie_polygons,
               aes(x = x, y = y, group = Category, fill = Category),
               colour = "white", linewidth = 0.35) +
  geom_text(data = label_data,
            aes(x = x, y = y, label = Label),
            size = 5.2, family = "Arial", colour = "black") +
  scale_fill_manual(values = category_colors, drop = FALSE,
                    breaks = category_labels, name = NULL) +
  coord_equal(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05),
              expand = FALSE, clip = "off") +
  labs(tag = "A") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             keywidth = grid::unit(0.55, "cm"),
                             keyheight = grid::unit(0.35, "cm"))) +
  theme_void(base_size = 12, base_family = "Arial") +
  theme(
    plot.tag = element_text(size = 28, face = "bold", hjust = 0,
                            vjust = 1, family = "Arial"),
    plot.tag.position = c(0.03, 0.98),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 13, family = "Arial"),
    legend.margin = margin(0, 0, 4, 0),
    plot.margin = margin(12, 10, 8, 10)
  )

ggsave(file.path(base_dir, "Figure6A_occurrence_distribution.pdf"),
       occurrence_plot, width = 5.5, height = 5.5, units = "in",
       device = grDevices::cairo_pdf)
ggsave(file.path(base_dir, "Figure6A_occurrence_distribution.png"),
       occurrence_plot, width = 5.5, height = 5.5, units = "in", dpi = 600)

message("Figure6A plot completed: ", base_dir)
