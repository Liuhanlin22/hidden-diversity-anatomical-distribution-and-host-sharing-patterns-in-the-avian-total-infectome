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

input_file <- file.path(script_dir, "Transmission_Stacked_Bar_Data_final.csv")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

summary_data <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "Transmission_Level",
  "k__Bacteria",
  "k__Eukaryota",
  "k__Viruses"
)
missing_columns <- setdiff(required_columns, names(summary_data))
if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

level_order <- c(
  "Within-species",
  "Cross-species",
  "Cross-genus",
  "Cross-family",
  "Cross-order"
)
stack_order <- c("k__Viruses", "k__Bacteria", "k__Eukaryota")
class_labels <- c(
  "k__Viruses" = "Virus",
  "k__Bacteria" = "Bacteria",
  "k__Eukaryota" = "Eukaryota"
)
class_colors <- c(
  "k__Viruses" = "#DFA3A7",
  "k__Bacteria" = "#C9D9E8",
  "k__Eukaryota" = "#C4DDBF"
)

summary_data$Transmission_Level <- factor(
  summary_data$Transmission_Level,
  levels = level_order
)
summary_data <- summary_data[!is.na(summary_data$Transmission_Level), , drop = FALSE]
summary_data <- summary_data[order(summary_data$Transmission_Level), , drop = FALSE]

if (sum(as.matrix(summary_data[, stack_order])) != 460) {
  stop("Figure6C input must represent 460 unique species.")
}

plot_data <- do.call(rbind, lapply(seq_len(nrow(summary_data)), function(i) {
  data.frame(
    Transmission_Level = as.character(summary_data$Transmission_Level[i]),
    Class = stack_order,
    Count = as.numeric(unlist(summary_data[i, stack_order], use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}))
plot_data$Count[is.na(plot_data$Count)] <- 0
plot_data$x <- match(plot_data$Transmission_Level, level_order)
plot_data$Class <- factor(plot_data$Class, levels = stack_order)
plot_data$ymin <- 0
plot_data$ymax <- 0

for (x_value in sort(unique(plot_data$x))) {
  idx <- which(plot_data$x == x_value)
  idx <- idx[order(match(as.character(plot_data$Class[idx]), stack_order))]
  cumulative <- c(0, cumsum(plot_data$Count[idx]))
  plot_data$ymin[idx] <- cumulative[-length(cumulative)]
  plot_data$ymax[idx] <- cumulative[-1]
}

label_y <- function(value) {
  ifelse(value %% 100 == 0, as.character(value), "")
}

transmission_plot <- ggplot(plot_data) +
  geom_rect(
    aes(
      xmin = x - 0.31,
      xmax = x + 0.31,
      ymin = ymin,
      ymax = ymax,
      fill = Class
    ),
    colour = NA
  ) +
  scale_fill_manual(
    values = class_colors,
    breaks = stack_order,
    labels = unname(class_labels[stack_order]),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq_along(level_order),
    labels = level_order,
    limits = c(0.45, length(level_order) + 0.55),
    minor_breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 400, by = 100),
    minor_breaks = seq(0, 400, by = 50),
    labels = label_y,
    limits = c(0, 400),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL, fill = NULL, tag = "C") +
  theme_gray(base_size = 8, base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "#EBEBEB", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_line(colour = "white", linewidth = 0.55),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "white", linewidth = 0.55),
    panel.grid.minor.y = element_line(colour = "white", linewidth = 0.35),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "#B3B3B3", linewidth = 0.3),
    axis.ticks.length = grid::unit(1.5, "pt"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "#4D4D4D", size = 7),
    axis.text.y = element_text(colour = "#4D4D4D", size = 7),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.key.width = grid::unit(0.30, "cm"),
    legend.key.height = grid::unit(0.27, "cm"),
    legend.text = element_text(colour = "black", size = 8.2),
    legend.spacing.y = grid::unit(1.5, "pt"),
    plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 1),
    plot.tag.position = c(0.02, 0.99),
    plot.margin = margin(5, 5, 7, 5)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(script_dir, "Figure6C_transmission_stacked_bar.pdf"),
  plot = transmission_plot,
  width = 3.84,
  height = 2.4,
  units = "in",
  device = grDevices::cairo_pdf
)
ggsave(
  filename = file.path(script_dir, "Figure6C_transmission_stacked_bar.png"),
  plot = transmission_plot,
  width = 3.84,
  height = 2.4,
  units = "in",
  dpi = 600
)

message("Figure6C plot completed: ", script_dir)

