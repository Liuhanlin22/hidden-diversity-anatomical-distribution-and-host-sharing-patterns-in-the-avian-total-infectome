#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
setwd(script_dir)

depth_files <- sort(list.files(pattern = "\\.depth\\.txt$", full.names = TRUE))
annotation_file <- "coverage_annotation.tsv"

if (length(depth_files) == 0) {
  stop("No .depth.txt files found in the working directory.")
}
if (!file.exists(annotation_file)) {
  stop("Missing annotation file: ", annotation_file)
}

annotation <- read.delim(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("library", "virus", "reads", "RPM")
if (!all(required_columns %in% colnames(annotation))) {
  stop("Annotation file must contain: ", paste(required_columns, collapse = ", "))
}

read_depth <- function(path) {
  x <- read.delim(
    path,
    header = FALSE,
    col.names = c("contig", "position", "depth"),
    colClasses = c("character", "numeric", "numeric"),
    comment.char = "",
    stringsAsFactors = FALSE
  )
  if (nrow(x) == 0) stop("Empty depth file: ", basename(path))
  x
}

get_annotation <- function(library_id) {
  row <- annotation[annotation$library == library_id, , drop = FALSE]
  if (nrow(row) != 1) stop("Expected one annotation row for: ", library_id)
  row
}

draw_coverage <- function(path) {
  x <- read_depth(path)
  library_id <- sub("\\.depth\\.txt$", "", basename(path))
  info <- get_annotation(library_id)
  ymax <- max(x$depth, na.rm = TRUE)
  xmax <- max(x$position, na.rm = TRUE)
  ylim <- c(0, if (ymax > 0) ymax * 1.08 else 1)

  plot(x$position, x$depth,
       type = "n", ylim = ylim,
       xlab = "Position (bp)", ylab = "Depth",
       main = info$virus,
       font.main = 4, cex.main = 1.15,
       font.lab = 2, cex.lab = 1.0,
       font.axis = 1, cex.axis = 0.85)
  abline(h = pretty(c(0, ymax)), col = "grey85", lty = 2)
  polygon(c(x$position[1], x$position, tail(x$position, 1)),
          c(0, x$depth, 0), col = "#648AB4", border = NA)
  lines(x$position, x$depth, col = "#648AB4", lwd = 0.4)
  text(
    x = xmax * 0.98,
    y = ylim[2] * 0.96,
    labels = paste0("reads = ", info$reads, "\nRPM = ", formatC(info$RPM, format = "f", digits = 2)),
    adj = c(1, 1), cex = 1.05
  )
  box()
}

for (path in depth_files) {
  library_id <- sub("\\.depth\\.txt$", "", basename(path))
  pdf(file.path(script_dir, paste0(library_id, "_coverage.pdf")),
      width = 8, height = 4, useDingbats = FALSE)
  par(mar = c(4.5, 4.8, 3.5, 1.2), mgp = c(2.7, 0.8, 0), las = 1)
  draw_coverage(path)
  dev.off()
}

n <- length(depth_files)
ncol <- min(4, n)
nrow <- ceiling(n / ncol)
combined_file <- file.path(script_dir, paste0("combined_", nrow, "x", ncol, "_coverage.pdf"))

pdf(combined_file, width = 6 * ncol, height = 4 * nrow, useDingbats = FALSE)
par(mfrow = c(nrow, ncol), mar = c(4.5, 4.8, 3.5, 1.2), mgp = c(2.7, 0.8, 0), las = 1)
for (path in depth_files) draw_coverage(path)
if (n < nrow * ncol) {
  for (i in seq_len(nrow * ncol - n)) plot.new()
}
dev.off()

message("Completed coverage plots: ", n)
message("Combined output: ", combined_file)
