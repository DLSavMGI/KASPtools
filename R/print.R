#' @export
print.kasp_run <- function(x, ...) {
  cat("<kasp_run>\n")
  cat("  Protocol: ", x$protocol_name, "\n", sep = "")
  cat("  Source:   ", x$source_type, "\n", sep = "")
  cat("  File:     ", x$source_file, "\n", sep = "")
  cat("  Wells:    ", nrow(x$plate_map), "\n", sep = "")
  cat("  Dyes:     ", paste(x$dyes, collapse = ", "), "\n", sep = "")
  cat("  Acq.:     ", paste(range(x$data$acquisition, na.rm = TRUE), collapse = "-"), "\n", sep = "")
  invisible(x)
}

#' @export
print.kasp_result <- function(x, ...) {
  cat("<kasp_result>\n")
  cat("  Protocol:      ", x$protocol_name, "\n", sep = "")
  cat("  Normalization: ", x$normalization, "\n", sep = "")
  cat("  Plateau:       ", paste(range(x$plateau), collapse = "-"), "\n", sep = "")
  cat("  Endpoint:      ", x$endpoint_acquisition, "\n", sep = "")
  cat("  Samples:       ", nrow(x$samples), "\n", sep = "")
  tab <- table(x$samples$Genotype, useNA = "ifany")
  if (length(tab) > 0) print(tab)
  invisible(x)
}
