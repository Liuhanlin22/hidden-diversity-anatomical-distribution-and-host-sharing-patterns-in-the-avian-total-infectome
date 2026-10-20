suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

input_file <- file.path(script_dir, "Figure1BC_input_counts.xlsx")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

read_count_sheet <- function(path, sheet) {
  dat <- read_excel(path, sheet = sheet, col_names = FALSE)
  dat <- dat[, 1:2]
  colnames(dat) <- c("category", "n")
  dat$category <- trimws(as.character(dat$category))
  dat$n <- as.numeric(dat$n)
  dat[!is.na(dat$category) & !is.na(dat$n), ]
}

order_data <- read_count_sheet(input_file, 2)
sample_data <- read_count_sheet(input_file, 1)
order_data$category[order_data$category == "pooling"] <- "Mixed species"
sample_data$category[sample_data$category == "Fecal samples"] <- "Fecal sample"
order_data$category <- factor(order_data$category, levels = order_data$category)
sample_data$category <- factor(sample_data$category, levels = sample_data$category)

order_colors <- colorRampPalette(c(
  "#641011", "#8B1D25", "#A9232B", "#C93435", "#DF5140",
  "#E8755A", "#EE987F", "#F2B4A0", "#F5D6CA", "#F7E5DE"
))(nrow(order_data))

sample_colors <- colorRampPalette(c(
  "#18396F", "#1E5795", "#2E72A8", "#478AB8", "#65A0C3",
  "#89B2CD", "#A8C3D9", "#C6D6E5", "#DCE6EF", "#E9EFF5"
))(nrow(sample_data))

base_theme <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    axis.title = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 8.5, colour = "black", angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8.5, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(colour = "#D9D9D9", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
    plot.margin = margin(6, 8, 8, 8)
  )

make_bar_plot <- function(dat, colors, tag) {
  ggplot(dat, aes(x = category, y = n, fill = category)) +
    geom_col(width = 0.72, colour = NA) +
    geom_text(aes(label = n, y = n * 1.18), size = 3.0, family = "Arial", vjust = 0) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_y_log10(
      limits = c(1, 1200),
      breaks = c(1, 10, 100, 1000),
      labels = c("1", "10", "100", "1000"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = NULL, y = "Number of libraries", tag = tag) +
    base_theme
}

plot_b <- make_bar_plot(order_data, order_colors, "B")
plot_c <- make_bar_plot(sample_data, sample_colors, "C")
figure_bc <- plot_b + plot_c + plot_layout(ncol = 2, widths = c(1, 1))

ggsave(file.path(script_dir, "Figure1BC_host_library_composition.pdf"), figure_bc,
       width = 12.3, height = 3.35, units = "in", device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(script_dir, "Figure1BC_host_library_composition.png"), figure_bc,
       width = 12.3, height = 3.35, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(file.path(script_dir, "Figure1B_order_composition.pdf"), plot_b,
       width = 6.1, height = 3.35, units = "in", device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(script_dir, "Figure1C_sample_type_composition.pdf"), plot_c,
       width = 6.1, height = 3.35, units = "in", device = cairo_pdf, limitsize = FALSE)
