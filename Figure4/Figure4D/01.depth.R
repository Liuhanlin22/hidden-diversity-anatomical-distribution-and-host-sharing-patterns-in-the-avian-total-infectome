#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
setwd(script_dir)

bam_files <- list.files(pattern = "\\.sorted\\.bam$", full.names = TRUE)
if (length(bam_files) == 0) {
  message("No BAM files found; existing depth files will be used.")
  quit(save = "no", status = 0)
}

if (Sys.which("samtools") == "") {
  stop("samtools is required to convert BAM files to depth files.")
}

for (bam_file in bam_files) {
  depth_file <- sub("\\.sorted\\.bam$", ".depth.txt", bam_file)
  if (file.exists(depth_file)) {
    message("Skipping existing file: ", basename(depth_file))
    next
  }

  status <- system2("samtools", c("depth", bam_file), stdout = depth_file)
  if (status != 0) stop("samtools depth failed: ", basename(bam_file))
  message("Created: ", basename(depth_file))
}
