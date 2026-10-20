required_packages <- "igraph"
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(igraph))

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

node_type <- V(network)$node_type
node_color <- c(
  "Host (Order)" = "#B8DDF5",
  "k__Viruses" = "#D81B60",
  "k__Bacteria" = "#C5A6D1",
  "k__Eukaryota" = "#94C983"
)
node_size <- ifelse(node_type == "Host (Order)", 4.2, 2.0)
label_size <- ifelse(node_type == "Host (Order)", 0.54, 0.40)
label_font <- ifelse(node_type == "Host (Order)", 2, 1)
node_fill <- unname(node_color[node_type])

draw_network <- function() {
  par(mar = c(0, 0, 0, 0), xpd = NA, family = "sans")
  plot(
    network,
    layout = network_layout,
    rescale = TRUE,
    asp = 0,
    margin = 0,
    vertex.color = node_fill,
    vertex.size = node_size,
    vertex.frame.color = "white",
    vertex.frame.width = 0.35,
    vertex.label = V(network)$name,
    vertex.label.color = "black",
    vertex.label.cex = label_size,
    vertex.label.family = "sans",
    vertex.label.font = label_font,
    edge.color = grDevices::adjustcolor("grey45", alpha.f = 0.28),
    edge.width = 0.25,
    edge.arrow.size = 0
  )
  legend(
    "topleft",
    legend = c("Host (order)", "Virus", "Bacteria", "Eukaryota"),
    pch = 21,
    pt.bg = c("#B8DDF5", "#D81B60", "#C5A6D1", "#94C983"),
    col = "#D9D9D9",
    pt.cex = c(1.55, 0.95, 0.95, 0.95),
    pt.lwd = 0.45,
    bty = "o",
    box.lwd = 0.45,
    box.col = "#777777",
    bg = grDevices::adjustcolor("white", alpha.f = 0.94),
    cex = 0.92,
    inset = c(0.015, 0.085),
    x.intersp = 0.9,
    y.intersp = 1.05
  )
  plot_usr <- par("usr")
  text(
    plot_usr[1] + 0.008 * diff(plot_usr[1:2]),
    plot_usr[4] - 0.008 * diff(plot_usr[3:4]),
    labels = "D",
    adj = c(0, 1),
    cex = 1.25,
    font = 2,
    family = "sans"
  )
}

grDevices::cairo_pdf(
  file.path(script_dir, "Figure6D_host_pathogen_network.pdf"),
  width = 8.4,
  height = 5.1,
  family = "sans",
  bg = "transparent"
)
draw_network()
grDevices::dev.off()

grDevices::png(
  file.path(script_dir, "Figure6D_host_pathogen_network.png"),
  width = 5040,
  height = 3060,
  res = 600,
  bg = "transparent"
)
draw_network()
grDevices::dev.off()

message("Figure6D network completed: ", script_dir)
