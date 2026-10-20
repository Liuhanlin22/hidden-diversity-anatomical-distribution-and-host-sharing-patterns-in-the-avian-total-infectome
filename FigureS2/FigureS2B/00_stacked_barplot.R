# FigureS2B: Family composition by host order
# R implementation of 00_stacked_barplot.py, aligned to the supplied
# Supplement_Figure2_qc_info.pdf reference panel B.

options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("This script requires the 'readxl' package. Install it with install.packages('readxl').")
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    normalizePath(dirname(sub("^--file=", "", file_arg[[1]])),
                   winslash = "\\", mustWork = FALSE)
  } else {
    normalizePath(getwd(), winslash = "\\", mustWork = FALSE)
  }
}

script_dir <- get_script_dir()
input_file <- file.path(script_dir, "00_qc.xlsx")
pdf_file <- file.path(script_dir, "Family_Composition_by_Order_R.pdf")
png_file <- file.path(script_dir, "Family_Composition_by_Order_R.png")
csv_file <- file.path(script_dir, "Family_Composition_by_Order_R_data.csv")

if (!file.exists(input_file)) {
  stop("Input workbook not found: ", input_file)
}

qc <- as.data.frame(
  readxl::read_excel(input_file, sheet = 1, .name_repair = "minimal"),
  check.names = FALSE
)
names(qc) <- trimws(names(qc))

required_columns <- c("Order", "Family")
missing_columns <- setdiff(required_columns, names(qc))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

qc[["Order"]] <- trimws(as.character(qc[["Order"]]))
qc[["Family"]] <- trimws(as.character(qc[["Family"]]))
qc <- qc[nzchar(qc[["Order"]]) & nzchar(qc[["Family"]]) &
           !is.na(qc[["Order"]]) & !is.na(qc[["Family"]]), , drop = FALSE]

# The supplied reference panel excludes the pooling category.  Remove it
# from both dimensions so it cannot appear as an extra bar or legend entry.
qc <- qc[
  tolower(qc[["Order"]]) != "pooling" &
    tolower(qc[["Family"]]) != "pooling",
  ,
  drop = FALSE
]

order_levels <- c(
  "Accipitriformes",
  "Anseriformes",
  "Apterygiformes",
  "Charadriiformes",
  "Columbiformes",
  "Falconiformes",
  "Galliformes",
  "Passeriformes",
  "Pelecaniformes",
  "Procellariiformes",
  "Psittaciformes",
  "Sphenisciformes",
  "Strigiformes",
  "Suliformes"
)

family_levels <- c(
  "Acanthizidae",
  "Accipitridae",
  "Anatidae",
  "Apterygidae",
  "Ardeidae",
  "Cacatuidae",
  "Callaeidae",
  "Cardinalidae",
  "Charadriidae",
  "Columbidae",
  "Falconidae",
  "Fringillidae",
  "Hirundinidae",
  "Icteridae",
  "Laridae",
  "Locustellidae",
  "Meliphagidae",
  "Mimidae",
  "Mohouidae",
  "Muscicapidae",
  "Oceanitidae",
  "Parulidae",
  "Passerellidae",
  "Petroicidae",
  "Phasianidae",
  "Procellariidae",
  "Prunellidae",
  "Psittaculidae",
  "Recurvirostridae",
  "Rhipiduridae",
  "Scolopacidae",
  "Spheniscidae",
  "Strigidae",
  "Strigopidae",
  "Sturnidae",
  "Sulidae",
  "Turdidae",
  "Vireonidae",
  "Zosteropidae"
)

unexpected_orders <- setdiff(unique(qc[["Order"]]), order_levels)
unexpected_families <- setdiff(unique(qc[["Family"]]), family_levels)
if (length(unexpected_orders) > 0L) {
  stop("Unexpected Order value(s): ", paste(unexpected_orders, collapse = "; "))
}
if (length(unexpected_families) > 0L) {
  stop("Unexpected Family value(s): ", paste(unexpected_families, collapse = "; "))
}

# Count records within each Order x Family combination, then normalize each
# Order column to 100%, matching pandas crosstab(normalize = "index").
count_matrix <- table(
  factor(qc[["Family"]], levels = family_levels),
  factor(qc[["Order"]], levels = order_levels)
)
percentage_matrix <- sweep(count_matrix, 2, colSums(count_matrix), "/") * 100

# Exact colors sampled from the supplied reference PDF's Family legend.
family_colors <- c(
  "Acanthizidae" = "#993404",
  "Accipitridae" = "#CC4C02",
  "Anatidae" = "#EC7014",
  "Apterygidae" = "#FEC44F",
  "Ardeidae" = "#FFF7BC",
  "Cacatuidae" = "#CB181D",
  "Callaeidae" = "#EF3B2C",
  "Cardinalidae" = "#FB6A4A",
  "Charadriidae" = "#FC9272",
  "Columbidae" = "#FEE0D2",
  "Falconidae" = "#810F7C",
  "Fringillidae" = "#88419D",
  "Hirundinidae" = "#8C6BB1",
  "Icteridae" = "#8C96C6",
  "Laridae" = "#9EBCDA",
  "Locustellidae" = "#E0ECF4",
  "Meliphagidae" = "#004529",
  "Mimidae" = "#006837",
  "Mohouidae" = "#238443",
  "Muscicapidae" = "#41AB5D",
  "Oceanitidae" = "#78C679",
  "Parulidae" = "#ADDE8E",
  "Passerellidae" = "#D9F0A3",
  "Petroicidae" = "#084081",
  "Phasianidae" = "#0868AC",
  "Procellariidae" = "#2B8CBE",
  "Prunellidae" = "#4EB3D3",
  "Psittaculidae" = "#74A9CF",
  "Recurvirostridae" = "#A6BDDB",
  "Rhipiduridae" = "#7A0177",
  "Scolopacidae" = "#AE017E",
  "Spheniscidae" = "#DD3497",
  "Strigidae" = "#F768A1",
  "Strigopidae" = "#FA9FB5",
  "Sturnidae" = "#FDE0DD",
  "Sulidae" = "#8DD3C7",
  "Turdidae" = "#E0E0E0",
  "Vireonidae" = "#BABABA",
  "Zosteropidae" = "#878787"
)

if (!identical(names(family_colors), family_levels)) {
  stop("Family color mapping is incomplete or out of order.")
}

# Export the exact table used for plotting.
plot_data <- data.frame(
  Order = rep(order_levels, each = length(family_levels)),
  Family = rep(family_levels, times = length(order_levels)),
  Percentage = as.vector(percentage_matrix),
  stringsAsFactors = FALSE
)
plot_data <- plot_data[plot_data$Percentage > 0, , drop = FALSE]
plot_data$Color <- unname(family_colors[plot_data$Family])
write.csv(plot_data, csv_file, row.names = FALSE, na = "")

draw_family_plot <- function() {
  n_order <- length(order_levels)
  x <- seq_len(n_order)
  bar_width <- 0.70

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(0.72, 0.28))

  # Main 100% stacked bar panel.
  par(mar = c(6.4, 4.9, 1.1, 0.2), family = "sans", xpd = NA)
  plot(NA, xlim = c(0.5, n_order + 0.5), ylim = c(0, 100),
       axes = FALSE, xlab = "", ylab = "", bty = "n")

  bottom <- rep(0, n_order)
  for (family_name in family_levels) {
    heights <- as.numeric(percentage_matrix[family_name, ])
    rect(x - bar_width / 2, bottom,
         x + bar_width / 2, bottom + heights,
         col = family_colors[[family_name]], border = NA)
    bottom <- bottom + heights
  }

  axis(2, at = seq(0, 100, by = 20), las = 1, lwd = 0.8, cex.axis = 0.9)
  axis(1, at = x, labels = FALSE, lwd = 0.8, tck = -0.015)
  text(x, y = -8.0, labels = order_levels, srt = 45,
       adj = c(1, 1), cex = 0.82, font = 1, xpd = NA)
  mtext("Percentage (%)", side = 2, line = 3.0, font = 2, cex = 1.0)
  mtext("Order", side = 1, line = 5.2, font = 2, cex = 1.0)
  mtext("B", side = 3, line = 0.05, adj = 0, font = 2, cex = 1.35)
  segments(0.5, 0, n_order + 0.5, 0, col = "#333333", lwd = 0.8)

  # Two-column Family legend, with italic family names as in the reference.
  par(mar = c(0, 0, 0, 0), family = "sans", xpd = NA)
  plot.new()
  legend("topleft", legend = family_levels,
         fill = unname(family_colors[family_levels]), border = NA,
         title = "Family", bty = "n", ncol = 2,
         cex = 0.78, text.font = 3, title.font = 2,
         y.intersp = 1.05, x.intersp = 0.45, inset = c(0.01, 0.01))
}

grDevices::cairo_pdf(pdf_file, width = 20, height = 7, bg = "transparent")
draw_family_plot()
grDevices::dev.off()

grDevices::png(png_file, width = 20, height = 7, units = "in", res = 300,
               bg = "transparent", type = "cairo")
draw_family_plot()
grDevices::dev.off()

message("FigureS2B R plot completed: ", pdf_file)
message("Filtered plotting table: ", csv_file)

