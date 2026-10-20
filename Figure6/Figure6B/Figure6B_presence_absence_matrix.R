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
input_file <- file.path(base_dir, "virus_widespread_species_details_final.csv")
if (!file.exists(input_file)) stop("Missing input file: ", input_file)

country_order <- c(
  "Antarctica", "Australia", "Brazil", "China", "Germany",
  "New Zealand", "the United States"
)

species_order <- c(
  "Gentoo penguin virus",
  "Campylobacter jejuni",
  "Hydatigera sp. Avi1",
  "Trichomonas sp. Avi2"
)

present_color <- "#405B92"
absent_color <- "#F0F0F0"
segment_present_color <- "#405B92"

details <- read.csv(input_file, check.names = FALSE,
                    stringsAsFactors = FALSE)
required_columns <- c("kindom", "family", "species_name", country_order)
missing_columns <- setdiff(required_columns, names(details))
if (length(missing_columns) > 0) {
  stop("Input is missing: ", paste(missing_columns, collapse = ", "))
}

details$species_name <- trimws(as.character(details$species_name))
for (country in country_order) {
  details[[country]] <- suppressWarnings(as.numeric(details[[country]]))
}
if (anyNA(details[, country_order])) {
  stop("Country abundance columns must contain numeric values.")
}
if (any(details[, country_order] < 0, na.rm = TRUE)) {
  stop("Country abundance values must be non-negative.")
}
if (anyDuplicated(details$species_name)) {
  stop("species_name must be unique in the input file.")
}

species_order <- c(species_order,
                   setdiff(details$species_name, species_order))
missing_species <- setdiff(species_order, details$species_name)
if (length(missing_species) > 0) {
  species_order <- setdiff(species_order, missing_species)
}
details <- details[match(species_order, details$species_name), , drop = FALSE]
details$species_name <- factor(details$species_name, levels = rev(species_order))

presence_matrix <- details[, country_order, drop = FALSE] > 0
presence_output <- details
presence_output[, country_order] <- lapply(
  presence_matrix,
  as.integer
)
write.csv(
  presence_output,
  file.path(base_dir, "virus_widespread_species_presence_absence_final.csv"),
  row.names = FALSE
)

row_positions <- setNames(seq_along(species_order), rev(species_order))
long_data <- do.call(
  rbind,
  lapply(seq_len(nrow(details)), function(i) {
    data.frame(
      species_name = as.character(details$species_name[i]),
      country = country_order,
      x = seq_along(country_order),
      y = row_positions[as.character(details$species_name[i])],
      present = as.logical(unlist(presence_matrix[i, country_order],
                                  use.names = FALSE)),
      stringsAsFactors = FALSE
    )
  })
)
long_data$point_color <- ifelse(long_data$present,
                                present_color, absent_color)

segment_data <- do.call(
  rbind,
  Filter(Negate(is.null), lapply(seq_along(species_order), function(i) {
    row_values <- as.logical(unlist(presence_matrix[i, country_order],
                                    use.names = FALSE))
    present_positions <- which(row_values)
    if (length(present_positions) < 2) return(NULL)
    segment_start <- seq.int(min(present_positions),
                             max(present_positions) - 1)
    data.frame(
      species_name = species_order[i],
      x = segment_start,
      xend = segment_start + 1,
      y = row_positions[species_order[i]],
      yend = row_positions[species_order[i]],
      segment_color = segment_present_color,
      stringsAsFactors = FALSE
    )
  }))
)

species_labels <- data.frame(
  species_name = species_order,
  x = 0.82,
  y = unname(row_positions[species_order]),
  stringsAsFactors = FALSE
)
country_labels <- data.frame(
  country = country_order,
  x = seq_along(country_order),
  y = length(species_order) + 0.55,
  stringsAsFactors = FALSE
)

presence_plot <- ggplot() +
  geom_point(data = long_data,
             aes(x = x, y = y, colour = point_color),
             size = 6.7, shape = 16) +
  geom_segment(data = segment_data,
               aes(x = x, xend = xend, y = y, yend = yend,
                   colour = segment_color),
               linewidth = 1.35, lineend = "round") +
  geom_text(data = species_labels,
            aes(x = x, y = y, label = species_name),
            hjust = 1, fontface = "italic", family = "Arial",
            size = 4.5, colour = "black") +
  geom_text(data = country_labels,
            aes(x = x, y = y, label = country),
            angle = 45, hjust = 0, vjust = 0,
            family = "Arial", size = 4.5, colour = "black") +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0.65, length(country_order) + 0.35),
                     breaks = seq_along(country_order),
                     labels = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.35, length(species_order) + 1.45),
                     breaks = NULL, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 12, base_family = "Arial") +
  theme(
    plot.margin = margin(18, 12, 8, 160),
    legend.position = "none"
  )

# Compact layout for the four-row panel.
figure_height <- 3.5

ggsave(file.path(base_dir, "Figure6B_presence_absence_matrix.pdf"),
       presence_plot, width = 9.6, height = figure_height, units = "in",
       device = grDevices::cairo_pdf)
ggsave(file.path(base_dir, "Figure6B_presence_absence_matrix.png"),
       presence_plot, width = 9.6, height = figure_height, units = "in", dpi = 600)

message("Figure6B plot completed: ", base_dir)
