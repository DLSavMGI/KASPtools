#' Export KASP analysis results
#'
#' Writes genotype calls, controls, model settings and optional technical
#' coordinates to an XLSX workbook. The interactive allelic-discrimination
#' plot can also be saved as HTML.
#'
#' @param result A `kasp_result` object returned by `kasp_cartesian()`.
#' @param xlsx_output Path to the output XLSX workbook.
#' @param html_output Optional path to save the interactive Plotly graph.
#' @param include_coordinates Include the full technical coordinate table.
#' @param overwrite Overwrite an existing XLSX file.
#'
#' @return Invisibly returns a list of created paths.
#' @export
export_kasp_results <- function(
  result,
  xlsx_output,
  html_output = NULL,
  include_coordinates = FALSE,
  overwrite = TRUE
) {
  if (!inherits(result, "kasp_result") && !inherits(result, "dt_kasp_result")) {
    stop("Expected an object returned by kasp_cartesian().", call. = FALSE)
  }
  if (missing(xlsx_output) || length(xlsx_output) != 1L || is.na(xlsx_output) || xlsx_output == "") {
    stop("xlsx_output must be supplied.", call. = FALSE)
  }

  samples_export <- result$samples |>
    dplyr::transmute(
      well,
      well_position,
      sample_name,
      Genotype,
      Confidence = dplyr::if_else(
        Confidence_Category == "Review",
        "Review",
        paste0(Confidence_Category, " (", Call_Confidence_Score, "%)")
      ),
      Note = Review_Note
    )

  controls_export <- result$controls |>
    dplyr::transmute(
      well,
      well_position,
      sample_name,
      task,
      Genotype,
      Confidence = "Control",
      Note = ""
    )

  ntc_export <- result$ntc |>
    dplyr::transmute(
      well,
      well_position,
      sample_name,
      task,
      Genotype = "NTC",
      Confidence = "NTC",
      Note = ""
    )

  all_calls <- dplyr::bind_rows(samples_export, controls_export, ntc_export)
  all_calls <- .add_plate_order(all_calls, "well_position")
  all_calls <- all_calls[order(
    is.na(all_calls$.plate_column),
    all_calls$.plate_column,
    all_calls$.plate_row,
    all_calls$well
  ), , drop = FALSE]
  all_calls$.plate_column <- NULL
  all_calls$.plate_row <- NULL
  names(all_calls)[names(all_calls) == "well"] <- "Well"
  names(all_calls)[names(all_calls) == "well_position"] <- "Well Position"
  names(all_calls)[names(all_calls) == "sample_name"] <- "Sample Name"
  names(all_calls)[names(all_calls) == "task"] <- "Task"

  controls_all <- dplyr::bind_rows(controls_export, ntc_export)
  controls_all <- .add_plate_order(controls_all, "well_position")
  controls_all <- controls_all[order(
    is.na(controls_all$.plate_column),
    controls_all$.plate_column,
    controls_all$.plate_row,
    controls_all$well
  ), , drop = FALSE]
  controls_all$.plate_column <- NULL
  controls_all$.plate_row <- NULL
  names(controls_all)[names(controls_all) == "well"] <- "Well"
  names(controls_all)[names(controls_all) == "well_position"] <- "Well Position"
  names(controls_all)[names(controls_all) == "sample_name"] <- "Sample Name"
  names(controls_all)[names(controls_all) == "task"] <- "Task"

  genotype_hetero <- unname(result$call_to_genotype["Heterozygote"])
  settings_rows <- c(
    Protocol = result$protocol_name,
    Source = result$source_type,
    Source_file = result$source_file,
    Normalization = result$normalization,
    Allele_1 = result$allele_1_name,
    Allele_2 = result$allele_2_name,
    Heterozygote = genotype_hetero,
    Allele_1_dye = result$allele_1_dye,
    Allele_2_dye = result$allele_2_dye,
    ROX_dye = result$rox_dye,
    Endpoint_acquisition = result$endpoint_acquisition,
    Plateau_acquisitions = paste0(min(result$plateau), "-", max(result$plateau)),
    Positive_signal_reference = result$positive_signal_reference
  )
  validation_rows <- unlist(result$settings, recursive = TRUE, use.names = TRUE)
  settings_export <- data.frame(
    Parameter = c(names(settings_rows), names(validation_rows)),
    Value = as.character(c(settings_rows, validation_rows)),
    stringsAsFactors = FALSE
  )

  wb <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9EAF7",
    border = "Bottom",
    borderColour = "#5B7C99"
  )
  review_style <- openxlsx::createStyle(fgFill = "#FCE4D6")
  high_style <- openxlsx::createStyle(fgFill = "#E2F0D9")
  moderate_style <- openxlsx::createStyle(fgFill = "#FFF2CC")

  .write_table <- function(sheet, data, table_name) {
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeDataTable(wb, sheet = sheet, x = data, tableName = table_name)
    openxlsx::addStyle(
      wb,
      sheet = sheet,
      style = header_style,
      rows = 1,
      cols = seq_len(ncol(data)),
      gridExpand = TRUE
    )
    openxlsx::freezePane(wb, sheet = sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet = sheet, cols = seq_len(ncol(data)), widths = "auto")
  }

  .write_table("genotype_calls", all_calls, "genotype_calls")
  .write_table("controls_and_NTC", controls_all, "controls_and_NTC")
  .write_table("cluster_model", result$cluster_model, "cluster_model")
  .write_table("settings", settings_export, "settings")

  confidence_column <- match("Confidence", names(all_calls))
  if (!is.na(confidence_column)) {
    review_rows <- which(all_calls$Confidence == "Review") + 1L
    high_rows <- which(grepl("^High", all_calls$Confidence)) + 1L
    moderate_rows <- which(grepl("^Moderate", all_calls$Confidence)) + 1L
    if (length(review_rows) > 0) {
      openxlsx::addStyle(wb, "genotype_calls", review_style, review_rows, seq_len(ncol(all_calls)), gridExpand = TRUE, stack = TRUE)
    }
    if (length(high_rows) > 0) {
      openxlsx::addStyle(wb, "genotype_calls", high_style, high_rows, seq_len(ncol(all_calls)), gridExpand = TRUE, stack = TRUE)
    }
    if (length(moderate_rows) > 0) {
      openxlsx::addStyle(wb, "genotype_calls", moderate_style, moderate_rows, seq_len(ncol(all_calls)), gridExpand = TRUE, stack = TRUE)
    }
  }

  if (isTRUE(include_coordinates)) {
    coordinates <- result$coordinates
    .write_table("coordinates", coordinates, "coordinates")
  }

  openxlsx::saveWorkbook(wb, file = xlsx_output, overwrite = overwrite)

  if (!is.null(html_output)) {
    tryCatch(
      htmlwidgets::saveWidget(result$plot, file = html_output, selfcontained = TRUE),
      error = function(e) {
        warning("Could not save self-contained HTML; saving with dependencies directory.", call. = FALSE)
        htmlwidgets::saveWidget(result$plot, file = html_output, selfcontained = FALSE)
      }
    )
  }

  invisible(list(xlsx = xlsx_output, html = html_output))
}
