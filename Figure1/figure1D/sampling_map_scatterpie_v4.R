suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scatterpie)
  library(ggnewscale)
  library(maps)
  library(grid)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

input_file <- file.path(script_dir, "data_v4_input.csv")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

data <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    Origin_raw = Origin,
    Origin = case_when(
      Origin %in% c("America", "United States", "USA") ~ "USA",
      Origin %in% c("United Kingdom", "UK") ~ "UK",
      Origin %in% c("Czechoslovakia", "Czech Republic") ~ "Czech Republic",
      Origin %in% c("Taiwan", "China") ~ "China",
      TRUE ~ Origin
    ),
    Form = case_when(
      Form %in% c("Wild", "Domestic") ~ Form,
      TRUE ~ as.character(Form)
    ),
    Latitude = if_else(Origin == "Antarctica", -64.1, as.numeric(Latitude)),
    Longitude = if_else(Origin == "Antarctica", -60.2, as.numeric(Longitude))
  )

country_counts <- data %>%
  count(Origin, Form, name = "Count") %>%
  pivot_wider(names_from = Form, values_from = Count, values_fill = 0)

if (!"Wild" %in% names(country_counts)) country_counts$Wild <- 0
if (!"Domestic" %in% names(country_counts)) country_counts$Domestic <- 0

country_coordinates <- data %>%
  group_by(Origin) %>%
  summarise(
    Latitude = first(Latitude),
    Longitude = first(Longitude),
    .groups = "drop"
  )

wider_data <- country_counts %>%
  left_join(country_coordinates, by = "Origin") %>%
  mutate(
    Total = Wild + Domestic,
    r = pmax(0.8, pmin(12, log10(Total + 1) * 3))
  ) %>%
  select(Origin, Latitude, Longitude, Wild, Domestic, Total, r)

if (nrow(wider_data) != n_distinct(wider_data$Origin)) {
  stop("Country-level aggregation failed.")
}

map_colors <- c("#fed976", "#fd8d3c", "#e31a1c")
type_colors <- c("Wild" = "#9699CA", "Domestic" = "#CBDDC7")

world <- map_data("world") %>%
  mutate(region = case_when(
    region == "Taiwan" ~ "China",
    region == "United States" ~ "USA",
    TRUE ~ region
  ))

world_data <- world %>%
  left_join(wider_data, by = c("region" = "Origin"))

base_map <- ggplot() +
  geom_polygon(
    data = world_data,
    aes(x = long, y = lat, group = group, fill = Total),
    color = "grey35",
    linewidth = 0.2
  ) +
  scale_fill_gradientn(
    colours = map_colors,
    values = c(0, 0.5, 1),
    limits = c(0, 300),
    breaks = c(0, 100, 200, 300),
    labels = c("0", "100", "200", "300"),
    na.value = "grey90",
    name = "Total samples"
  ) +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-82, 83), expand = FALSE) +
  labs(tag = "D") +
  theme_void(base_family = "Arial") +
  theme(
    plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
    plot.tag.position = c(0.01, 0.99),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.height = unit(0.55, "cm"),
    panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.35),
    plot.margin = margin(6, 8, 6, 8)
  )

max_total <- max(wider_data$Total, na.rm = TRUE)
if (max_total != 203) stop("Unexpected maximum sample count.")
legend_values <- c(0, 6, 35, 203)
legend_radii <- log10(legend_values + 1) * 3

figure_d <- base_map +
  ggnewscale::new_scale_fill() +
  geom_scatterpie(
    data = wider_data,
    aes(x = Longitude, y = Latitude, r = r),
    cols = c("Wild", "Domestic"),
    color = "black",
    linewidth = 0.3,
    alpha = 1
  ) +
  scale_fill_manual(values = type_colors, name = "Type") +
  geom_scatterpie_legend(
    r = legend_radii,
    x = -157,
    y = -69,
    breaks = legend_radii,
    labeller = function(x) round(10^(x / 3) - 1, 0),
    size = 3,
    n = length(legend_values)
  ) +
  annotate("text", x = -157, y = -53, label = "Sample size", size = 4, fontface = "bold", hjust = 0)

ggsave(
  file.path(script_dir, "sampling_map_wild_domestic_v4.pdf"),
  figure_d,
  width = 12,
  height = 6,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  file.path(script_dir, "sampling_map_wild_domestic_v4.png"),
  figure_d,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
