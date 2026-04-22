#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript convert_seqlogo_rds_to_json.R <input_seqlogo.rds> <output_seqlogo.json>")
}

input_path <- args[[1]]
output_path <- args[[2]]

seqlogos <- readRDS(input_path)

if (!is.list(seqlogos) || is.null(names(seqlogos))) {
  stop("Expected seqlogo.rds to be a named list of motif probability matrices.")
}

bases <- c("A", "C", "G", "T")

seqlogos_json <- lapply(names(seqlogos), function(motif_name) {
  m <- seqlogos[[motif_name]]

  if (!is.matrix(m)) {
    stop(sprintf("Motif '%s' is not a matrix.", motif_name))
  }

  if (!all(bases %in% rownames(m))) {
    stop(sprintf(
      "Motif '%s' is missing one or more required rownames: A, C, G, T.",
      motif_name
    ))
  }

  # R object is bases x positions. Python-friendly JSON should be
  # positions x bases, as a plain list of dicts:
  # [{"A": ..., "C": ..., "G": ..., "T": ...}, ...]
  df <- as.data.frame(t(m[bases, , drop = FALSE]))
  colnames(df) <- bases

  # Ensure numeric values and no row names leak into JSON.
  df[] <- lapply(df, as.numeric)
  rownames(df) <- NULL

  df
})

names(seqlogos_json) <- names(seqlogos)

write_json(
  seqlogos_json,
  path = output_path,
  pretty = TRUE,
  auto_unbox = TRUE,
  dataframe = "rows",
  digits = NA
)

message(sprintf("Wrote %d motif logos to %s", length(seqlogos_json), output_path))

