.kasp_expand_program <- function(program_node) {
  step_nodes <- xml2::xml_find_all(program_node, "./step")
  steps <- purrr::map_dfr(step_nodes, function(step_node) {
    temperature_node <- xml2::xml_find_first(step_node, "./temperature")
    loop_node <- xml2::xml_find_first(step_node, "./loop")
    step_type <- dplyr::case_when(
      !inherits(temperature_node, "xml_missing") ~ "temperature",
      !inherits(loop_node, "xml_missing") ~ "loop",
      TRUE ~ "other"
    )
    tibble::tibble(
      program_step = suppressWarnings(as.integer(.kasp_xml_text(step_node, "./nr"))),
      step_type = step_type,
      temperature_C = suppressWarnings(as.numeric(.kasp_xml_text(step_node, "./temperature/temperature"))),
      duration_s = suppressWarnings(as.numeric(.kasp_xml_text(step_node, "./temperature/duration"))),
      measure_mode = .kasp_xml_text(step_node, "./temperature/measure"),
      goto_step = suppressWarnings(as.integer(.kasp_xml_text(step_node, "./loop/goto"))),
      repeat_n = suppressWarnings(as.integer(.kasp_xml_text(step_node, "./loop/repeat")))
    )
  }) |>
    dplyr::arrange(program_step)

  loop_definitions <- steps |>
    dplyr::filter(step_type == "loop")

  loop_iterations <- stats::setNames(
    rep(1L, nrow(loop_definitions)),
    as.character(loop_definitions$program_step)
  )

  acquisition_events <- list()
  step_index <- 1L
  acquisition_index <- 0L
  safety_counter <- 0L

  while (step_index <= nrow(steps)) {
    safety_counter <- safety_counter + 1L
    if (safety_counter > 100000L) {
      stop("Could not expand the thermal program: too many transitions.", call. = FALSE)
    }

    current <- steps[step_index, , drop = FALSE]

    if (current$step_type == "temperature") {
      active_loops <- loop_definitions |>
        dplyr::filter(goto_step <= current$program_step, program_step > current$program_step)

      if (nrow(active_loops) > 0) {
        active_loops <- active_loops |>
          dplyr::mutate(loop_width = program_step - goto_step) |>
          dplyr::arrange(loop_width)
        inner_loop <- active_loops[1, , drop = FALSE]
        loop_step <- inner_loop$program_step
        loop_repeat <- inner_loop$repeat_n
        cycle_in_block <- loop_iterations[[as.character(loop_step)]]
      } else {
        loop_step <- NA_integer_
        loop_repeat <- NA_integer_
        cycle_in_block <- 1L
      }

      if (!is.na(current$measure_mode) && current$measure_mode != "") {
        acquisition_index <- acquisition_index + 1L
        acquisition_events[[length(acquisition_events) + 1L]] <- tibble::tibble(
          acquisition = acquisition_index,
          program_step = current$program_step,
          temperature_C = current$temperature_C,
          duration_s = current$duration_s,
          loop_step = loop_step,
          cycle_in_block = cycle_in_block,
          repeat_in_block = loop_repeat
        )
      }
      step_index <- step_index + 1L

    } else if (current$step_type == "loop") {
      loop_key <- as.character(current$program_step)
      current_iteration <- loop_iterations[[loop_key]]

      if (current_iteration < current$repeat_n) {
        loop_iterations[[loop_key]] <- current_iteration + 1L

        nested_loops <- loop_definitions |>
          dplyr::filter(
            program_step < current$program_step,
            goto_step >= current$goto_step
          )
        if (nrow(nested_loops) > 0) {
          for (nested_step in nested_loops$program_step) {
            loop_iterations[[as.character(nested_step)]] <- 1L
          }
        }

        goto_index <- match(current$goto_step, steps$program_step)
        if (is.na(goto_index)) {
          stop("Thermal program refers to a missing goto step: ", current$goto_step, call. = FALSE)
        }
        step_index <- goto_index
      } else {
        step_index <- step_index + 1L
      }

    } else {
      step_index <- step_index + 1L
    }
  }

  program <- dplyr::bind_rows(acquisition_events)
  if (nrow(program) == 0) {
    stop("No fluorescence acquisition steps were found in the thermal program.", call. = FALSE)
  }

  program |>
    dplyr::mutate(
      cycle_block = match(program_step, unique(program_step)),
      cycle_name = dplyr::if_else(
        !is.na(repeat_in_block),
        paste0(
          "Block ", cycle_block, " | ", temperature_C, " C | ",
          cycle_in_block, "/", repeat_in_block
        ),
        paste0("Block ", cycle_block, " | ", temperature_C, " C")
      )
    )
}

#' Import DT96 RDML/XML fluorescence data
#'
#' Imports all fluorescence channels and all acquisition points from a DT96
#' RDML/XML export. Separate assay blocks can be selected automatically by
#' gaps between occupied plate columns or explicitly with `columns`.
#'
#' @param file Path to an RDML/XML file.
#' @param protocol_to_keep Integer index of the automatically detected assay block.
#' @param protocol_name User-facing assay/protocol name stored in the returned object.
#' @param columns Optional integer vector of plate columns to select explicitly.
#' @param output_rds Optional path for saving the imported object as RDS.
#'
#' @return A `kasp_run` object.
#' @export
import_kasp_xml <- function(
  file,
  protocol_to_keep = 1L,
  protocol_name,
  columns = NULL,
  output_rds = NULL
) {
  if (missing(protocol_name) || length(protocol_name) != 1L || is.na(protocol_name) || trimws(protocol_name) == "") {
    stop("protocol_name must be supplied.", call. = FALSE)
  }
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)

  xml <- xml2::read_xml(file)
  xml2::xml_ns_strip(xml)

  target_nodes <- xml2::xml_find_all(xml, "./target")
  target_map <- purrr::map_dfr(target_nodes, function(x) {
    tibble::tibble(
      target_id = xml2::xml_attr(x, "id"),
      dye = .kasp_xml_attr(x, "./dyeId", "id")
    )
  }) |>
    dplyr::mutate(
      target_id = trimws(as.character(target_id)),
      dye = toupper(trimws(as.character(dye)))
    )

  program_nodes <- xml2::xml_find_all(xml, "./thermalCyclingConditions")
  if (length(program_nodes) == 0) {
    stop("No thermalCyclingConditions were found in the RDML file.", call. = FALSE)
  }
  if (length(program_nodes) > 1) {
    warning("Multiple thermalCyclingConditions were found; the first program is used.", call. = FALSE)
  }
  program <- .kasp_expand_program(program_nodes[[1]])

  run_nodes <- xml2::xml_find_all(xml, "./experiment/run")
  if (length(run_nodes) == 0) stop("No run elements were found in the RDML file.", call. = FALSE)

  raw_data <- purrr::map_dfr(run_nodes, function(run_node) {
    plate_columns <- suppressWarnings(as.integer(.kasp_xml_text(run_node, "./pcrFormat/columns")))
    if (is.na(plate_columns)) stop("Could not determine the number of plate columns.", call. = FALSE)

    react_nodes <- xml2::xml_find_all(run_node, "./react")
    purrr::map_dfr(react_nodes, function(react_node) {
      react_id <- xml2::xml_attr(react_node, "id")
      well_index <- suppressWarnings(as.integer(react_id))
      if (is.na(well_index)) {
        stop("Non-numeric react id is not supported by this importer: ", react_id, call. = FALSE)
      }

      row_index <- ((well_index - 1L) %/% plate_columns) + 1L
      column_index <- ((well_index - 1L) %% plate_columns) + 1L
      if (row_index > length(LETTERS)) stop("Plate row exceeds supported range.", call. = FALSE)
      well_position <- paste0(LETTERS[row_index], column_index)
      sample_name <- .kasp_xml_text(react_node, "./sample")

      data_nodes <- xml2::xml_find_all(react_node, "./data")
      purrr::map_dfr(data_nodes, function(data_node) {
        current_target_id <- .kasp_xml_attr(data_node, "./tar", "id")
        dye <- target_map |>
          dplyr::filter(target_id == current_target_id) |>
          dplyr::pull(dye)
        if (length(dye) == 0 || is.na(dye[1]) || dye[1] == "") {
          dye <- toupper(current_target_id)
        } else {
          dye <- dye[1]
        }

        adp_nodes <- xml2::xml_find_all(data_node, "./adp")
        if (length(adp_nodes) == 0) return(tibble::tibble())
        fluorescence <- suppressWarnings(as.numeric(
          xml2::xml_text(xml2::xml_find_first(adp_nodes, "./fluor"))
        ))

        tibble::tibble(
          well = well_index,
          well_position = well_position,
          row_index = row_index,
          column_index = column_index,
          sample_name = sample_name,
          dye = dye,
          acquisition = seq_along(adp_nodes),
          fluorescence = fluorescence
        )
      })
    })
  })

  if (nrow(raw_data) == 0) stop("No fluorescence data were found in the XML file.", call. = FALSE)

  sample_conflicts <- raw_data |>
    dplyr::distinct(well, sample_name) |>
    dplyr::count(well, name = "n_sample_names") |>
    dplyr::filter(n_sample_names > 1)
  if (nrow(sample_conflicts) > 0) {
    stop("Different sample names were found for the same well across runs.", call. = FALSE)
  }

  plate_map_all <- raw_data |>
    dplyr::distinct(well, well_position, row_index, column_index, sample_name)
  occupied_columns <- sort(unique(plate_map_all$column_index))
  protocol_columns <- tibble::tibble(column_index = occupied_columns) |>
    dplyr::mutate(
      new_protocol = c(TRUE, diff(column_index) > 1),
      protocol_index = cumsum(new_protocol)
    )
  protocol_summary <- plate_map_all |>
    dplyr::left_join(
      protocol_columns |> dplyr::select(column_index, protocol_index),
      by = "column_index"
    ) |>
    dplyr::group_by(protocol_index) |>
    dplyr::summarise(
      columns = paste(sort(unique(column_index)), collapse = ", "),
      wells = dplyr::n_distinct(well),
      .groups = "drop"
    )

  if (is.null(columns)) {
    protocol_to_keep <- as.integer(protocol_to_keep)
    if (length(protocol_to_keep) != 1L || is.na(protocol_to_keep) ||
        !protocol_to_keep %in% protocol_summary$protocol_index) {
      stop(
        "protocol_to_keep was not found. Available blocks: ",
        paste(protocol_summary$protocol_index, collapse = ", "),
        call. = FALSE
      )
    }
    selected_columns <- protocol_columns |>
      dplyr::filter(protocol_index == protocol_to_keep) |>
      dplyr::pull(column_index)
  } else {
    selected_columns <- sort(unique(as.integer(columns)))
    selected_columns <- selected_columns[is.finite(selected_columns)]
    if (length(selected_columns) == 0) stop("columns did not contain valid plate columns.", call. = FALSE)
    missing_columns <- setdiff(selected_columns, occupied_columns)
    if (length(missing_columns) > 0) {
      stop("Requested plate columns are not occupied: ", paste(missing_columns, collapse = ", "), call. = FALSE)
    }
    protocol_to_keep <- NA_integer_
  }

  kasp_data <- raw_data |>
    dplyr::filter(column_index %in% selected_columns) |>
    dplyr::left_join(program, by = "acquisition") |>
    dplyr::transmute(
      well,
      well_position,
      sample_name,
      dye = toupper(trimws(dye)),
      acquisition,
      fluorescence,
      temperature_C,
      program_step,
      cycle_block,
      cycle_in_block,
      cycle_name
    ) |>
    dplyr::arrange(well, dye, acquisition)

  duplicates <- kasp_data |>
    dplyr::count(well, dye, acquisition) |>
    dplyr::filter(n > 1)
  if (nrow(duplicates) > 0) {
    stop("Duplicate well x dye x acquisition records were found.", call. = FALSE)
  }

  plate_map <- kasp_data |>
    dplyr::distinct(well, well_position, sample_name) |>
    dplyr::mutate(task = "UNKNOWN", omit = FALSE) |>
    dplyr::arrange(well)

  result <- list(
    protocol_name = trimws(protocol_name),
    protocol_index = protocol_to_keep,
    selected_columns = selected_columns,
    source_file = basename(file),
    source_type = "DT96_RDML",
    dyes = sort(unique(kasp_data$dye)),
    protocol_summary = protocol_summary,
    plate_map = plate_map,
    data = kasp_data
  )
  class(result) <- c("kasp_run", "dt_kasp", "list")

  if (!is.null(output_rds)) saveRDS(result, file = output_rds)
  result
}
