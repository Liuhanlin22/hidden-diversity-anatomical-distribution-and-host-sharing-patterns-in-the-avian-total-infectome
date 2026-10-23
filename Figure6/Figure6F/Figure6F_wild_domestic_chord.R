required_packages <- "circlize"
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(circlize))

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  getwd()
}

input_file <- file.path(script_dir, "virus_table_merged_by_Domestic_or_Wild_updated.csv")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_columns <- c("species_name", "kindom", "Domestic", "Wild")
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

data$Domestic <- suppressWarnings(as.numeric(data$Domestic))
data$Wild <- suppressWarnings(as.numeric(data$Wild))
data$Domestic[is.na(data$Domestic)] <- 0
data$Wild[is.na(data$Wild)] <- 0
data$species_name <- trimws(data$species_name)
data$kindom <- trimws(data$kindom)

data <- data[data$Domestic > 0 & data$Wild > 0, , drop = FALSE]
if (nrow(data) == 0) {
  stop("No species has positive abundance in both Domestic and Wild.")
}
if (nrow(data) != 22) {
  stop("The final Figure6F input must contain 22 taxa detected in both groups.")
}
if (anyDuplicated(data$species_name)) {
  stop("species_name must be unique in the Domestic/Wild input table.")
}
if (any(grepl("^Alphainfluenzavirus influenzae \\(", data$species_name)) ||
    any(data$species_name %in% c(
      "Streptococcus suis (strain avian group)",
      "Streptococcus suis (strain pig group)"
    ))) {
  stop("Influenza and Streptococcus suis labels must be merged at species level.")
}

write.table(
  data.frame(
    species_name = data$species_name,
    kindom = data$kindom,
    Domestic_RPM_sum = data$Domestic,
    Wild_RPM_sum = data$Wild,
    stringsAsFactors = FALSE
  ),
  file.path(script_dir, "Figure6F_input_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

kingdom_order <- c("k__Bacteria", "k__Eukaryota", "k__Viruses")
data$kindom <- factor(data$kindom, levels = kingdom_order)
data <- data[order(data$kindom, seq_len(nrow(data))), , drop = FALSE]
data$kindom <- as.character(data$kindom)

species <- data$species_name
species_kindom <- setNames(data$kindom, species)

link_data <- rbind(
  data.frame(
    from = "Domestic",
    to = species,
    value = log10(data$Domestic + 1),
    stringsAsFactors = FALSE
  ),
  data.frame(
    from = "Wild",
    to = species,
    value = log10(data$Wild + 1),
    stringsAsFactors = FALSE
  )
)

grid_colors <- c(
  Domestic = "#B9D5EB",
  Wild = "#FFFAEE",
  k__Viruses = "#E295A9",
  k__Bacteria = "#9B9CBF",
  k__Eukaryota = "#B5CFB0"
)
link_colors <- c(
  k__Viruses = "#F0CAD4",
  k__Bacteria = "#CDCDDF",
  k__Eukaryota = "#DAE7D7"
)
sector_colors <- c(
  Domestic = grid_colors[["Domestic"]],
  setNames(unname(grid_colors[species_kindom]), species),
  Wild = grid_colors[["Wild"]]
)
link_col <- vapply(link_data$to, function(x) {
  if (x %in% c("Domestic", "Wild")) {
    return("#D9D9D9")
  }
  grDevices::adjustcolor(link_colors[[species_kindom[[x]]]], alpha.f = 0.72)
}, character(1))

sector_order <- c(
  "Domestic",
  species[species_kindom == "k__Bacteria"],
  species[species_kindom == "k__Eukaryota"],
  species[species_kindom == "k__Viruses"],
  "Wild"
)
gap_after <- rep(1.5, length(sector_order))

grDevices::cairo_pdf(
  file.path(script_dir, "Figure6F_wild_domestic_chord.pdf"),
  width = 6.6,
  height = 6.6,
  family = "sans",
  bg = "transparent"
)
par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA, family = "sans")
circos.clear()
circos.par(
  # Rotate the sector sequence to match the final composite Figure 6F:
  # Domestic at upper-left, Wild at lower-left, and pathogens clockwise
  # through the right and lower arcs.
  start.degree = 170,
  canvas.xlim = c(-1.35, 1.35),
  canvas.ylim = c(-1.35, 1.35),
  gap.after = gap_after,
  track.margin = c(0.006, 0.006),
  cell.padding = c(0, 0, 0, 0),
  points.overflow.warning = FALSE
)

chordDiagram(
  link_data,
  order = sector_order,
  grid.col = sector_colors,
  col = link_col,
  directional = 0,
  transparency = 0,
  link.border = NA,
  link.lwd = 0.25,
  # Preserve the species order grouped by kingdom at both Domestic and Wild.
  # Sorting by abundance would interleave the three kingdoms.
  link.sort = FALSE,
  link.decreasing = FALSE,
  # The colored sector ring is drawn once below.  Do not also draw circlize's
  # default grid track, otherwise the figure contains two concentric rings.
  annotationTrack = NULL,
  # Match the narrow single outer sector ring in the reference Figure 6F.
  preAllocateTracks = list(track.height = 0.07),
  big.gap = 1.5,
  small.gap = 1.5
)

circos.trackPlotRegion(
  track.index = 1,
  ylim = c(0, 1),
  bg.border = NA,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    fill <- unname(sector_colors[sector_name])
    circos.rect(
      CELL_META$xlim[1], 0,
      CELL_META$xlim[2], 1,
      col = fill,
      border = "white",
      lwd = 0.35
    )
    if (sector_name == "Domestic") {
      circos.text(
        CELL_META$xcenter, 0.5, "Domestic",
        facing = "bending.inside",
        niceFacing = TRUE,
        col = "#2F5597",
        font = 2,
        cex = 0.85
      )
    } else if (sector_name == "Wild") {
      circos.text(
        CELL_META$xcenter, 0.5, "Wild",
        facing = "bending.inside",
        niceFacing = TRUE,
        col = "#D51F2A",
        font = 2,
        cex = 0.85
      )
    } else {
      circos.text(
        CELL_META$xcenter, 0.5, sector_name,
        facing = "clockwise",
        niceFacing = TRUE,
                  adj = c(0, 0.5),
                  cex = 0.50,
        col = "#222222"
      )
    }
  }
)

usr <- par("usr")
text(
  usr[1] + 0.01 * diff(usr[1:2]),
  usr[4] - 0.01 * diff(usr[3:4]),
  labels = "F",
  adj = c(0, 1),
  cex = 1.35,
  font = 2,
  family = "sans"
)
legend(
  "topleft",
  legend = c("Virus", "Bacteria", "Eukaryota"),
  fill = unname(grid_colors[c("k__Viruses", "k__Bacteria", "k__Eukaryota")]),
  border = NA,
  bty = "n",
  cex = 0.84,
  y.intersp = 0.9,
  x.intersp = 0.7,
  inset = c(0.055, 0.06)
)
dev.off()
circos.clear()

grDevices::png(
  file.path(script_dir, "Figure6F_wild_domestic_chord.png"),
  width = 3960,
  height = 3960,
  res = 600,
  bg = "transparent"
)
par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA, family = "sans")
circos.par(
  start.degree = 170,
  canvas.xlim = c(-1.35, 1.35),
  canvas.ylim = c(-1.35, 1.35),
  gap.after = gap_after,
  track.margin = c(0.006, 0.006),
  cell.padding = c(0, 0, 0, 0),
  points.overflow.warning = FALSE
)
chordDiagram(
  link_data,
  order = sector_order,
  grid.col = sector_colors,
  col = link_col,
  directional = 0,
  transparency = 0,
  link.border = NA,
  link.lwd = 0.25,
  link.sort = FALSE,
  link.decreasing = FALSE,
  # Keep a single custom colored ring, matching the reference figure.
  annotationTrack = NULL,
  preAllocateTracks = list(track.height = 0.07),
  big.gap = 1.5,
  small.gap = 1.5
)
circos.trackPlotRegion(
  track.index = 1,
  ylim = c(0, 1),
  bg.border = NA,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    fill <- unname(sector_colors[sector_name])
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                col = fill, border = "white", lwd = 0.35)
    if (sector_name == "Domestic") {
      circos.text(CELL_META$xcenter, 0.5, "Domestic", facing = "bending.inside",
                  niceFacing = TRUE, col = "#2F5597", font = 2, cex = 0.85)
    } else if (sector_name == "Wild") {
      circos.text(CELL_META$xcenter, 0.5, "Wild", facing = "bending.inside",
                  niceFacing = TRUE, col = "#D51F2A", font = 2, cex = 0.85)
    } else {
      circos.text(CELL_META$xcenter, 0.5, sector_name, facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5), cex = 0.50, col = "#222222")
    }
  }
)
usr <- par("usr")
text(usr[1] + 0.01 * diff(usr[1:2]), usr[4] - 0.01 * diff(usr[3:4]),
     labels = "F", adj = c(0, 1), cex = 1.35, font = 2, family = "sans")
legend("topleft", legend = c("Virus", "Bacteria", "Eukaryota"),
       fill = unname(grid_colors[c("k__Viruses", "k__Bacteria", "k__Eukaryota")]),
       border = NA, bty = "n", cex = 0.84, y.intersp = 0.9,
       x.intersp = 0.7, inset = c(0.055, 0.06))
dev.off()
circos.clear()

write.csv(
  data.frame(
    from = link_data$from,
    to = link_data$to,
    log10_RPM_plus_1 = link_data$value,
    link_class = vapply(link_data$to, function(x) {
      if (x %in% c("Domestic", "Wild")) "Host" else species_kindom[[x]]
    }, character(1)),
    stringsAsFactors = FALSE
  ),
  file.path(script_dir, "Figure6F_links_used.csv"),
  row.names = FALSE,
  na = ""
)

message("Figure6F chord diagram completed: ", script_dir)
