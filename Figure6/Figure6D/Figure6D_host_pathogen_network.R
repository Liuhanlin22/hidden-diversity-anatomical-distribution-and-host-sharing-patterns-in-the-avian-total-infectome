required_packages <- c("igraph", "ggplot2", "ggrepel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  getwd()
}

input_candidates <- c(
  "virus_table_filtered_updated.csv",
  "Figure6D_final_host_order_matrix.csv"
)
input_paths <- file.path(script_dir, input_candidates)
existing_inputs <- input_paths[file.exists(input_paths)]
if (length(existing_inputs) == 0) {
  stop("No Figure6D input matrix found in: ", script_dir)
}
input_file <- existing_inputs[1]

data <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(data) < 5 || !all(c("species_name", "kindom") %in% names(data))) {
  stop("The input must contain kindom, species_name, and host-order columns.")
}
if (nrow(data) != 60 || anyDuplicated(data$species_name)) {
  stop("The final Figure6D input must contain 60 unique species-level taxa.")
}
if (any(grepl("^Alphainfluenzavirus influenzae \\(", data$species_name)) ||
    any(data$species_name %in% c(
      "Streptococcus suis (strain avian group)",
      "Streptococcus suis (strain pig group)"
    ))) {
  stop("Influenza and Streptococcus suis labels must be merged at species level.")
}

host_cols <- setdiff(names(data)[5:ncol(data)], "pooling")
host_cols <- host_cols[nzchar(host_cols)]
host_cols <- host_cols[vapply(data[host_cols], function(values) {
  numeric_values <- suppressWarnings(as.numeric(values))
  any(!is.na(numeric_values) & numeric_values > 0)
}, logical(1))]
if (length(host_cols) == 0) {
  stop("No host-order columns found.")
}

edge_list <- do.call(rbind, lapply(seq_len(nrow(data)), function(i) {
  values <- suppressWarnings(as.numeric(data[i, host_cols]))
  present <- host_cols[!is.na(values) & values > 0]
  species <- trimws(as.character(data$species_name[i]))
  kindom <- trimws(as.character(data$kindom[i]))
  if (length(present) == 0 || is.na(species) || !nzchar(species)) {
    return(NULL)
  }
  data.frame(
    Microbe = rep(species, length(present)),
    Host = present,
    Kindom = rep(kindom, length(present)),
    stringsAsFactors = FALSE
  )
}))

if (is.null(edge_list) || nrow(edge_list) == 0) {
  stop("No positive host-pathogen associations were found.")
}

edge_list <- unique(edge_list)
write.csv(
  edge_list,
  file.path(script_dir, "Figure6D_network_edges_used.csv"),
  row.names = FALSE,
  na = ""
)
selection_audit <- data.frame(
  lib_id = data$lib_id,
  species_name = data$species_name,
  positive_host_orders = vapply(seq_len(nrow(data)), function(i) {
    sum(suppressWarnings(as.numeric(data[i, host_cols])) > 0, na.rm = TRUE)
  }, integer(1)),
  selection_reason = ">=2_host_orders",
  status = "retained",
  stringsAsFactors = FALSE
)
write.csv(
  selection_audit,
  file.path(script_dir, "Figure6D_node_selection_audit.csv"),
  row.names = FALSE,
  na = ""
)
graph_edges <- edge_list[c("Microbe", "Host")]
names(graph_edges) <- c("from", "to")

microbe_vertices <- unique(edge_list[c("Microbe", "Kindom")])
names(microbe_vertices) <- c("name", "node_type")
host_vertices <- data.frame(
  name = host_cols,
  node_type = "Host (Order)",
  stringsAsFactors = FALSE
)
vertices <- rbind(host_vertices, microbe_vertices)
vertices <- vertices[!duplicated(vertices$name), , drop = FALSE]

network <- graph_from_data_frame(graph_edges, directed = FALSE, vertices = vertices)
set.seed(42)
network_layout <- layout_with_fr(network, niter = 12000, start.temp = 0.35, grid = "nogrid")
host_index <- V(network)$node_type == "Host (Order)"
host_center <- colMeans(network_layout[host_index, , drop = FALSE])
network_layout[host_index, ] <- host_center + 0.78 * (network_layout[host_index, ] - host_center)
# Match the vertical orientation of the original Figure6D layout.
network_layout[, 2] <- -network_layout[, 2]
network_layout[, 1] <- -network_layout[, 1]
layout_rotation <- -pi / 6
rotation_matrix <- matrix(
  c(cos(layout_rotation), -sin(layout_rotation),
    sin(layout_rotation), cos(layout_rotation)),
  nrow = 2,
  byrow = TRUE
)
network_layout <- network_layout %*% rotation_matrix

node_type <- as.character(V(network)$node_type)
node_color <- c(
  "Host (Order)" = "#B8DDF5",
  "k__Viruses" = "#D81B60",
  "k__Bacteria" = "#C5A6D1",
  "k__Eukaryota" = "#94C983"
)
node_data <- data.frame(
  name = V(network)$name,
  node_type = node_type,
  x = network_layout[, 1],
  y = network_layout[, 2],
  stringsAsFactors = FALSE
)
node_data$label <- vapply(node_data$name, function(label) {
  paste(strwrap(label, width = 32), collapse = "\n")
}, character(1))
node_data$point_size <- ifelse(node_data$node_type == "Host (Order)", 3.5, 2.2)
node_data$label_size <- ifelse(node_data$node_type == "Host (Order)", 2.2, 1.55)
node_data$label_font <- ifelse(node_data$node_type == "Host (Order)", "bold", "plain")

edge_data <- as_data_frame(network, what = "edges")
edge_data$x <- node_data$x[match(edge_data$from, node_data$name)]
edge_data$y <- node_data$y[match(edge_data$from, node_data$name)]
edge_data$xend <- node_data$x[match(edge_data$to, node_data$name)]
edge_data$yend <- node_data$y[match(edge_data$to, node_data$name)]

network_plot <- ggplot() +
  geom_segment(
    data = edge_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = grDevices::adjustcolor("grey45", alpha.f = 0.28),
    linewidth = 0.25
  ) +
  geom_point(
    data = node_data,
    aes(x = x, y = y, fill = node_type, size = point_size),
    shape = 21,
    colour = "white",
    stroke = 0.35
  ) +
  geom_text_repel(
    data = node_data,
    aes(x = x, y = y, label = label, size = label_size, fontface = label_font),
    colour = "black",
    box.padding = grid::unit(0.30, "lines"),
    point.padding = grid::unit(0.12, "lines"),
    force = 5,
    force_pull = 0.005,
    max.time = 60,
    max.iter = 100000,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.colour = grDevices::adjustcolor("grey35", alpha.f = 0.55),
    segment.size = 0.16,
    seed = 42,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    name = NULL,
    values = node_color,
    breaks = c("Host (Order)", "k__Viruses", "k__Bacteria", "k__Eukaryota"),
    labels = c("Host (order)", "Virus", "Bacteria", "Eukaryota")
  ) +
  scale_size_identity(guide = "none") +
  coord_equal(clip = "off") +
  labs(tag = "D") +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.text = element_text(size = 7),
    legend.key.size = grid::unit(0.16, "in"),
    legend.spacing.x = grid::unit(0.08, "in"),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.005, 0.995),
    plot.margin = margin(8, 14, 8, 8),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  file.path(script_dir, "Figure6D_host_pathogen_network.pdf"),
  network_plot,
  width = 10.5,
  height = 6.3,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(script_dir, "Figure6D_host_pathogen_network.png"),
  network_plot,
  width = 10.5,
  height = 6.3,
  units = "in",
  dpi = 600,
  bg = "white"
)

message("Figure6D network completed: ", script_dir)
