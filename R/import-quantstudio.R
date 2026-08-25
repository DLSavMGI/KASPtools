#' Import QuantStudio KASP fluorescence data
#'
#' Reads the `Results` and `Multicomponent Data` sheets from a QuantStudio
#' XLSX export, selects one SNP assay, and converts the fluorescence data to
#' the same long-format object used by the DT96 importer.
#'
#' @param file Path to the QuantStudio XLSX export.
#' @param assay_to_keep Exact SNP assay name to select.
#' @param protocol_name User-facing assay name. Defaults to `assay_to_keep`.
#' @param dye_columns Candidate fluorescence columns to import when present.
#' @param include_omitted Whether wells marked Omit should be retained.
#' @param output_rds Optional path for saving the imported object as RDS.
#'
#' @return A `kasp_run` object.
#' @export
import_kasp_quantstudio <- function(
  file,
  assay_to_keep,
  protocol_name = assay_to_keep,
  dye_columns = c("FAM", "VIC", "HEX", "ROX"),
  include_omitted = FALSE,
  output_rds = NULL
) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  if (length(assay_to_keep) != 1L || is.na(assay_to_keep) || trimws(assay_to_keep) == "") {
    stop("assay_to_keep must be a single non-empty assay name.", call. = FALSE)
  }
  if (length(protocol_name) != 1L || is.na(protocol_name) || trimws(protocol_name) == "") {
    stop("protocol_name must be a single non-empty name.", call. = FALSE)
  }

  results_header_row <- .find_header_row(
    file_path = file,
    sheet_name = "Results",
    required_headers = c("Well", "Well Position", "Sample Name", "SNP Assay Name", "Task")
  )
  multicomponent_header_row <- .find_header_row(
    file_path = file,
    sheet_name = "Multicomponent Data",
    required_headers = c("Well", "Cycle")
  )

  results_raw <- .read_table_from_header(file, "Results", results_header_row)
  required_results_columns <- c("Well", "Well.Position", "Sample.Name", "SNP.Assay.Name", "Task")
  missing_results <- setdiff(required_results_columns, names(results_raw))
  if (length(missing_results) > 0) {
    stop("Results sheet is missing columns: ", paste(missing_results, collapse = ", "), call. = FALSE)
  }

  plate_map <- results_raw |>
    dplyr::filter(!is.na(Well)) |>
    dplyr::transmute(
      well = as.integer(.to_numeric(Well)),
      well_position = toupper(trimws(as.character(Well.Position))),
      sample_name = trimws(as.character(Sample.Name)),
      assay_name = trimws(as.character(SNP.Assay.Name)),
      task = toupper(trimws(as.character(Task))),
      omit = if ("Omit" %in% names(results_raw)) {
        tolower(trimws(as.character(Omit))) %in% c("true", "1", "yes")
      } else {
        FALSE
      }
    ) |>
    dplyr::mutate(sample_name = dplyr::na_if(sample_name, "")) |>
    dplyr::filter(assay_name == assay_to_keep)

  if (!isTRUE(include_omitted)) plate_map <- plate_map |> dplyr::filter(!omit)
  if (nrow(plate_map) == 0) {
    stop("Assay '", assay_to_keep, "' was not found or contains no usable wells.", call. = FALSE)
  }
  duplicate_wells <- plate_map |> dplyr::count(well) |> dplyr::filter(n > 1)
  if (nrow(duplicate_wells) > 0) {
    stop("Duplicate wells were found within the selected assay.", call. = FALSE)
  }

  fluorescence_raw <- .read_table_from_header(
    file,
    "Multicomponent Data",
    multicomponent_header_row
  )
  required_fluorescence_columns <- c("Well", "Cycle")
  missing_fluorescence <- setdiff(required_fluorescence_columns, names(fluorescence_raw))
  if (length(missing_fluorescence) > 0) {
    stop(
      "Multicomponent Data sheet is missing columns: ",
      paste(missing_fluorescence, collapse = ", "),
      call. = FALSE
    )
  }

  dye_columns <- unique(toupper(trimws(as.character(dye_columns))))
  available_dyes <- intersect(dye_columns, names(fluorescence_raw))
  if (length(available_dyes) < 2) {
    stop(
      "Fewer than two requested fluorescence channels were found. Present requested channels: ",
      paste(available_dyes, collapse = ", "),
      call. = FALSE
    )
  }

  fluorescence_wide <- fluorescence_raw |>
    dplyr::transmute(
      well = as.integer(.to_numeric(Well)),
      acquisition = as.integer(.to_numeric(Cycle)),
      dplyr::across(dplyr::all_of(available_dyes), .to_numeric)
    ) |>
    dplyr::filter(!is.na(well), !is.na(acquisition))

  fluorescence_long <- fluorescence_wide |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(available_dyes),
      names_to = "dye",
      values_to = "fluorescence"
    ) |>
    dplyr::filter(!is.na(fluorescence))

  data <- fluorescence_long |>
    dplyr::inner_join(
      plate_map |> dplyr::select(well, well_position, sample_name),
      by = "well"
    ) |>
    dplyr::transmute(
      well,
      well_position,
      sample_name,
      dye = toupper(dye),
      acquisition,
      fluorescence,
      temperature_C = NA_real_,
      program_step = NA_integer_,
      cycle_block = 1L,
      cycle_in_block = acquisition,
      cycle_name = paste0("Cycle ", acquisition)
    ) |>
    dplyr::arrange(well, dye, acquisition)

  result <- list(
    protocol_name = trimws(protocol_name),
    protocol_index = NA_integer_,
    selected_columns = sort(unique(suppressWarnings(as.integer(sub("^[A-Z]+", "", plate_map$well_position))))),
    source_file = basename(file),
    source_type = "QuantStudio_XLSX",
    dyes = sort(unique(data$dye)),
    protocol_summary = tibble::tibble(
      protocol_index = 1L,
      columns = paste(sort(unique(suppressWarnings(as.integer(sub("^[A-Z]+", "", plate_map$well_position))))), collapse = ", "),
      wells = nrow(plate_map)
    ),
    plate_map = plate_map |> dplyr::select(well, well_position, sample_name, task, omit),
    data = data
  )
  class(result) <- c("kasp_run", "quantstudio_kasp", "list")

  if (!is.null(output_rds)) saveRDS(result, file = output_rds)
  result
}
