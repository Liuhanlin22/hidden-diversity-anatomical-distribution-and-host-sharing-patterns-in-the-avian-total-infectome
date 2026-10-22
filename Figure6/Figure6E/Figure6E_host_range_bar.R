required_packages <- "ggplot2"
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(ggplot2))

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  getwd()
}

merged_file <- file.path(script_dir, "virus_table_merged_by_Order_updated.csv")

if (!file.exists(merged_file)) {
  stop("Input file not found: ", merged_file)
}

merged_data <- read.csv(
  merged_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_merged_columns <- c("species_name", "Accipitriformes", "Anseriformes",
                             "Apterygiformes", "Charadriiformes", "Columbiformes",
                             "Falconiformes", "Galliformes", "Passeriformes",
                             "Pelecaniformes", "Procellariiformes", "Psittaciformes",
                             "Sphenisciformes", "Strigiformes", "Suliformes")
missing_merged_columns <- setdiff(required_merged_columns, names(merged_data))
if (length(missing_merged_columns) > 0) {
  stop("Missing columns in merged input: ", paste(missing_merged_columns, collapse = ", "))
}

order_columns <- required_merged_columns[-1]
order_matrix <- as.data.frame(lapply(merged_data[order_columns], function(x) {
  value <- suppressWarnings(as.numeric(x))
  value[is.na(value)] <- 0
  value
}), check.names = FALSE)

plot_species <- merged_data$species_name
plot_species[plot_species %in% c(
  "Streptococcus suis (strain avian group)",
  "Streptococcus suis (strain pig group)"
)] <- "Streptococcus suis"
plot_species[grepl(
  "^Alphainfluenzavirus influenzae \\(",
  plot_species
)] <- "Alphainfluenzavirus influenzae"

presence_data <- data.frame(
  species_name = plot_species,
  order_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
presence_data <- aggregate(
  . ~ species_name,
  data = presence_data,
  FUN = sum
)

count_data <- data.frame(
  species_name = presence_data$species_name,
  cross_order_count = rowSums(presence_data[order_columns] > 0),
  stringsAsFactors = FALSE
)
count_data <- count_data[count_data$cross_order_count >= 3, , drop = FALSE]
count_data <- count_data[order(-count_data$cross_order_count, count_data$species_name), , drop = FALSE]
plot_data <- count_data
plot_data$x <- 1 + (seq_len(nrow(plot_data)) - 1) * 0.64
# Keep the compressed x-spacing, but narrow the bars independently.
# The plotting range is also compressed when x-spacing is reduced, so
# 0.56 would look almost unchanged in the exported figure. 0.42 gives
# an approximately 20% narrower visual bar than the previous output.
capsule_width <- 0.42
capsule_x_radius <- capsule_width / 2
capsule_y_radius <- 0.68

capsule_data <- do.call(rbind, lapply(seq_len(nrow(plot_data)), function(i) {
  x_center <- plot_data$x[i]
  y_max <- plot_data$cross_order_count[i]
  y_radius <- min(capsule_y_radius, y_max / 2)
  theta_bottom <- seq(pi, 2 * pi, length.out = 25)
  theta_top <- seq(0, pi, length.out = 25)
  bottom <- data.frame(
    x = x_center + capsule_x_radius * cos(theta_bottom),
    y = y_radius + y_radius * sin(theta_bottom)
  )
  right <- data.frame(
    x = c(x_center + capsule_x_radius, x_center + capsule_x_radius),
    y = c(y_radius, y_max - y_radius)
  )
  top <- data.frame(
    x = x_center + capsule_x_radius * cos(theta_top),
    y = y_max - y_radius + y_radius * sin(theta_top)
  )
  left <- data.frame(
    x = c(x_center - capsule_x_radius, x_center - capsule_x_radius),
    y = c(y_max - y_radius, y_radius)
  )
  data.frame(rbind(bottom, right, top, left), group = i)
}))

host_range_plot <- ggplot() +
  geom_hline(
    yintercept = seq(0, 8, by = 2),
    colour = "#E5E4E5",
    linewidth = 0.30,
    linetype = "solid"
  ) +
  geom_polygon(
    data = capsule_data,
    aes(x = x, y = y, group = group),
    fill = "#87C6CE",
    colour = NA
  ) +
  geom_text(
    data = plot_data,
    aes(x = x, y = cross_order_count, label = cross_order_count),
    vjust = -0.24,
    size = 2.35,
    family = "Arial",
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = plot_data$x,
    labels = plot_data$species_name,
    limits = c(0.60, max(plot_data$x) + 0.30),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 8, by = 2),
    limits = c(0, 8.5),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Host range (order)",
    tag = "E"
  ) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.35),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = grid::unit(1.5, "pt"),
    axis.text.x = element_text(
    angle = 55,
      hjust = 1,
      vjust = 1,
      size = 5.3,
      colour = "black"
    ),
    axis.text.y = element_text(size = 7, colour = "black"),
    axis.title.y = element_text(size = 8.5, colour = "black", margin = margin(r = 5)),
    plot.tag = element_text(size = 12, face = "bold", hjust = 0, vjust = 1),
    plot.tag.position = c(0.015, 0.99),
    plot.margin = margin(5, 5, 12, 8)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  file.path(script_dir, "Figure6E_host_range_bar.pdf"),
  host_range_plot,
  # Match the exported canvas width to the compressed x-spacing instead
  # of stretching the bars back across the original full width.
  width = 4.5,
  height = 2.4,
  units = "in",
  device = grDevices::cairo_pdf
)
ggsave(
  file.path(script_dir, "Figure6E_host_range_bar.png"),
  host_range_plot,
  width = 4.5,
  height = 2.4,
  units = "in",
  dpi = 600
)

message("Figure6E plot completed: ", script_dir)
