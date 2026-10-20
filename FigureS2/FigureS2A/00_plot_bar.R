# FigureS1A: QC up/down bar plot
# R implementation of 00_plot_bar.py.
# The output is written with an _R suffix so the existing Python PDF is preserved.

options(stringsAsFactors = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("This script requires the 'readxl' package. Install it with install.packages('readxl').")
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    normalizePath(dirname(sub("^--file=", "", file_arg[[1]])), winslash = "\\", mustWork = FALSE)
  } else {
    normalizePath(getwd(), winslash = "\\", mustWork = FALSE)
  }
}

script_dir <- get_script_dir()
input_file <- file.path(script_dir, "00_qc.xlsx")
pdf_file <- file.path(script_dir, "QC_up_down_barplot_ordered_R.pdf")
png_file <- file.path(script_dir, "QC_up_down_barplot_ordered_R.png")
csv_file <- file.path(script_dir, "QC_up_down_barplot_ordered_R_data.csv")

if (!file.exists(input_file)) {
  stop("Input workbook not found: ", input_file)
}

qc <- as.data.frame(
  readxl::read_excel(input_file, sheet = 1, .name_repair = "minimal"),
  check.names = FALSE
)
names(qc) <- trimws(names(qc))

required_columns <- c("lib_id", "Sample Type", "norRNA percentage", "rawdata_reads")
missing_columns <- setdiff(required_columns, names(qc))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

qc[["rawdata_reads"]] <- suppressWarnings(as.numeric(qc[["rawdata_reads"]]))
qc[["norRNA percentage"]] <- suppressWarnings(as.numeric(qc[["norRNA percentage"]]))
if (any(!is.finite(qc[["rawdata_reads"]])) || any(qc[["rawdata_reads"]] < 0)) {
  stop("'rawdata_reads' contains missing, non-numeric, or negative values.")
}
if (any(!is.finite(qc[["norRNA percentage"]]))) {
  stop("'norRNA percentage' contains missing or non-numeric values.")
}

qc[["log10_raw_reads"]] <- log10(qc[["rawdata_reads"]] + 1)

# This is the same biological/sample order used by the Python script.
sample_order_raw <- c(
  "Oral swab",
  "Cloacal swab",
  "Oral and cloacal swab",
  "Fecal samples",
  "Cecal content",
  "Brain",
  "Heart",
  "Spleen",
  "Kidney",
  "Lung",
  "Liver",
  "Intestine",
  "Heart and kidney",
  "Mixed lung, Liver, skin, heart muscle, breast muscle and Brain of a single chinck",
  "Lesion",
  "Blood",
  "Feathers"
)

# Short labels reproduce the labels in the supplied reference PDF while the
# original values remain available in the input-derived column.
display_labels <- c(
  "Heart and kidney" = "Mixed heart and kidney",
  "Mixed lung, Liver, skin, heart muscle, breast muscle and Brain of a single chinck" = "Mixed organ"
)
sample_order <- unname(ifelse(
  sample_order_raw %in% names(display_labels),
  display_labels[sample_order_raw],
  sample_order_raw
))

qc[["sample_type_raw"]] <- trimws(as.character(qc[["Sample Type"]]))
qc[["sample_type_plot"]] <- qc[["sample_type_raw"]]
recode_idx <- qc[["sample_type_plot"]] %in% names(display_labels)
qc[["sample_type_plot"]][recode_idx] <- unname(display_labels[qc[["sample_type_plot"]][recode_idx]])

unexpected <- setdiff(unique(qc[["sample_type_plot"]]), sample_order)
if (length(unexpected) > 0L) {
  stop("Unexpected Sample Type value(s): ", paste(unexpected, collapse = "; "))
}

qc[["sample_type_plot"]] <- factor(qc[["sample_type_plot"]], levels = sample_order, ordered = TRUE)
qc <- qc[order(as.integer(qc[["sample_type_plot"]]), -qc[["norRNA percentage"]], na.last = TRUE), , drop = FALSE]
row.names(qc) <- NULL

# Fixed colors reproduce the supplied reference PDF and remain stable across runs.
sample_colors <- c(
  "Oral swab" = "#A50F15",
  "Cloacal swab" = "#FB6A4A",
  "Oral and cloacal swab" = "#FC9272",
  "Fecal samples" = "#9ECAE1",
  "Cecal content" = "#DEEBF7",
  "Brain" = "#3F007D",
  "Heart" = "#6A51A3",
  "Spleen" = "#807DBA",
  "Kidney" = "#9E9AC8",
  "Lung" = "#BCBDDC",
  "Liver" = "#045A8D",
  "Intestine" = "#0570B0",
  "Mixed heart and kidney" = "#3690C0",
  "Mixed organ" = "#74A9CF",
  "Lesion" = "#238443",
  "Blood" = "#78C679",
  "Feathers" = "#D9F0A3"
)

qc[["color"]] <- unname(sample_colors[as.character(qc[["sample_type_plot"]])])
qc[["x_position"]] <- seq_len(nrow(qc))

# Save the exact plotting table for auditability and Illustrator/data checking.
write.csv(
  qc[, c("lib_id", "Sample Type", "sample_type_plot", "rawdata_reads",
         "norRNA percentage", "log10_raw_reads", "color", "x_position")],
  csv_file,
  row.names = FALSE,
  na = ""
)

draw_qc_plot <- function() {
  n <- nrow(qc)
  x <- qc[["x_position"]]
  bar_col <- qc[["color"]]

  top_max <- max(1, ceiling(max(qc[["norRNA percentage"]]) * 10) / 10)
  # Keep the lower limit close to the largest observed value.  Rounding this
  # limit up to 10 creates excessive white space; the reference panel ends
  # just below the largest bar (with 2, 4, 6 and 8 as the visible ticks).
  bottom_max <- max(1, max(qc[["log10_raw_reads"]]) * 1.02)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # The third layout cell is the one-column legend spanning the two plot rows.
  layout(matrix(c(1, 3, 2, 3), nrow = 2, byrow = TRUE),
         widths = c(0.79, 0.21), heights = c(1, 1))

  # Top panel: norRNA percentage.
  par(mar = c(0.1, 5.2, 1.15, 0.2), family = "sans", xpd = NA)
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0, top_max),
       axes = FALSE, xlab = "", ylab = "", bty = "n")
  rect(x - 0.5, 0, x + 0.5, qc[["norRNA percentage"]],
       col = bar_col, border = NA)
  axis(2, las = 1, lwd = 0.8, cex.axis = 0.9)
  mtext("norRNA percentage", side = 2, line = 3.0, font = 2, cex = 1.0)
  mtext("A", side = 3, line = 0.05, adj = 0, font = 2, cex = 1.35)
  # Draw the shared horizontal baseline only from the first to the last bar;
  # do not extend it through the empty half-bar margins.
  segments(1, 0, n, 0, col = "#333333", lwd = 0.8)

  # Bottom panel: log10(rawdata_reads), displayed below the shared zero line.
  par(mar = c(0.8, 5.2, 0.1, 0.2), family = "sans", xpd = NA)
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(-bottom_max, 0),
       axes = FALSE, xlab = "", ylab = "", bty = "n")
  rect(x - 0.5, 0, x + 0.5, -qc[["log10_raw_reads"]],
       col = bar_col, border = NA)
  bottom_ticks <- seq(2, bottom_max, by = 2)
  axis(2, at = c(0, -bottom_ticks), labels = c(0, bottom_ticks),
       las = 1, lwd = 0.8, cex.axis = 0.9)
  mtext("log10(Raw Reads)", side = 2, line = 3.0, font = 2, cex = 1.0)
  segments(1, 0, n, 0, col = "#333333", lwd = 0.8)

  # Legend panel, ordered exactly as the bars.
  par(mar = c(0, 0, 0, 0), family = "sans", xpd = NA)
  plot.new()
  legend("topleft", legend = sample_order,
         fill = unname(sample_colors[sample_order]), border = NA,
         title = "Sample Type", bty = "n", cex = 0.78,
         y.intersp = 1.05, x.intersp = 0.55, inset = c(0.01, 0.02))
}

grDevices::cairo_pdf(pdf_file, width = 8.4, height = 3.6, bg = "transparent")
draw_qc_plot()
grDevices::dev.off()

grDevices::png(png_file, width = 8.4, height = 3.6, units = "in", res = 300,
               bg = "transparent", type = "cairo")
draw_qc_plot()
grDevices::dev.off()

message("FigureS1A R plot completed: ", pdf_file)
message("Processed plotting table: ", csv_file)

