.kasp_xml_text <- function(node, xpath) {
  x <- xml2::xml_find_first(node, xpath)
  if (inherits(x, "xml_missing")) return(NA_character_)
  value <- trimws(xml2::xml_text(x))
  if (identical(value, "")) NA_character_ else value
}

.kasp_xml_attr <- function(node, xpath, attribute) {
  x <- xml2::xml_find_first(node, xpath)
  if (inherits(x, "xml_missing")) return(NA_character_)
  value <- xml2::xml_attr(x, attribute)
  if (is.na(value) || trimws(value) == "") NA_character_ else value
}

.to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
}

.find_header_row <- function(file_path, sheet_name, required_headers) {
  raw_sheet <- openxlsx::read.xlsx(
    xlsxFile = file_path,
    sheet = sheet_name,
    colNames = FALSE,
    detectDates = FALSE,
    skipEmptyRows = FALSE,
    skipEmptyCols = FALSE
  )
  header_row_found <- apply(raw_sheet, 1, function(x) {
    x <- trimws(as.character(x))
    all(required_headers %in% x)
  })
  possible_rows <- which(header_row_found)
  if (length(possible_rows) == 0) {
    stop(
      paste0(
        "Could not find the table header on sheet '", sheet_name,
        "'. Expected fields: ", paste(required_headers, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  possible_rows[1]
}

.read_table_from_header <- function(file_path, sheet_name, header_row) {
  dat <- openxlsx::read.xlsx(
    xlsxFile = file_path,
    sheet = sheet_name,
    startRow = header_row,
    colNames = TRUE,
    detectDates = FALSE,
    skipEmptyRows = TRUE,
    skipEmptyCols = FALSE,
    check.names = TRUE
  )
  names(dat) <- trimws(names(dat))
  dat
}

.normalize_wells <- function(x) {
  if (is.null(x)) return(character())
  x <- toupper(trimws(as.character(x)))
  x[!is.na(x) & x != ""]
}

.make_hover_text <- function(well_position, sample_name) {
  well_position <- trimws(as.character(well_position))
  sample_name <- trimws(as.character(sample_name))
  well_position[is.na(well_position) | well_position == ""] <- "not specified"
  sample_name[is.na(sample_name) | sample_name == ""] <- "not specified"
  paste0(
    "<b>Well:</b> ", well_position,
    "<br><b>Sample:</b> ", sample_name
  )
}

.clamp01 <- function(x) pmax(0, pmin(1, x))

.transform_to_control_basis <- function(x, y, ntc_x, ntc_y, basis_matrix) {
  shifted <- rbind(x - ntc_x, y - ntc_y)
  transformed <- solve(basis_matrix, shifted)
  tibble::tibble(
    A1_Component = as.numeric(transformed[1, ]),
    A2_Component = as.numeric(transformed[2, ])
  )
}

.make_axis_range <- function(values, padding_fraction = 0.03) {
  values <- values[is.finite(values)]
  if (length(values) == 0) stop("Could not determine plot axis range.", call. = FALSE)
  value_min <- min(values)
  value_max <- max(values)
  value_span <- value_max - value_min
  if (value_span < .Machine$double.eps) {
    padding <- max(abs(value_max) * padding_fraction, 0.01)
  } else {
    padding <- value_span * padding_fraction
  }
  c(value_min - padding, value_max + padding)
}

.add_plate_order <- function(data, well_position_column = "well_position") {
  well_position <- toupper(trimws(as.character(data[[well_position_column]])))
  valid_position <- grepl("^[A-Z]+[0-9]+$", well_position)
  plate_row <- ifelse(
    valid_position,
    sub("^([A-Z]+)([0-9]+)$", "\\1", well_position),
    NA_character_
  )
  plate_column <- ifelse(
    valid_position,
    suppressWarnings(as.integer(sub("^([A-Z]+)([0-9]+)$", "\\2", well_position))),
    NA_integer_
  )
  data$.plate_column <- plate_column
  data$.plate_row <- match(plate_row, LETTERS)
  data
}

.calculate_local_compactness <- function(data, k = 3L) {
  data$Local_kNN_Distance <- NA_real_
  data$Score_Local_Compactness <- NA_real_
  cluster_names <- unique(data$Empirical_Cluster[!is.na(data$Empirical_Cluster)])

  for (cluster_name in cluster_names) {
    cluster_index <- which(
      data$Empirical_Cluster == cluster_name &
        is.finite(data$A1_Component) &
        is.finite(data$A2_Component)
    )
    n_cluster <- length(cluster_index)
    if (n_cluster == 0) next
    if (n_cluster < 4) {
      data$Score_Local_Compactness[cluster_index] <- 0.80
      next
    }

    coordinates <- as.matrix(
      data[cluster_index, c("A1_Component", "A2_Component"), drop = FALSE]
    )
    distance_matrix <- as.matrix(stats::dist(coordinates))
    diag(distance_matrix) <- Inf
    k_current <- min(as.integer(k), n_cluster - 1L)
    local_distance <- apply(distance_matrix, 1, function(x) {
      mean(sort(x)[seq_len(k_current)])
    })
    median_distance <- stats::median(local_distance, na.rm = TRUE)
    distance_mad <- stats::mad(
      local_distance,
      center = median_distance,
      constant = 1,
      na.rm = TRUE
    )
    transition_scale <- max(3 * distance_mad, 0.50 * median_distance, 1e-6)
    local_score <- exp(-pmax(local_distance - median_distance, 0) / transition_scale)
    data$Local_kNN_Distance[cluster_index] <- local_distance
    data$Score_Local_Compactness[cluster_index] <- .clamp01(local_score)
  }

  data
}

.validate_kasp_run <- function(kasp) {
  if (!inherits(kasp, "kasp_run") && !inherits(kasp, "dt_kasp")) {
    stop("Expected an object created by a KASPtools import function.", call. = FALSE)
  }
  required <- c("protocol_name", "dyes", "plate_map", "data")
  missing <- setdiff(required, names(kasp))
  if (length(missing) > 0) {
    stop("Malformed KASP run object. Missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.control_task_to_call <- function() {
  c(
    "PC_ALLELE_1" = "Allele 1 / Allele 1",
    "PC_ALLELE_2" = "Allele 2 / Allele 2",
    "PC_ALLELE_BOTH" = "Heterozygote"
  )
}

.genotype_calls <- function() {
  c("Allele 1 / Allele 1", "Heterozygote", "Allele 2 / Allele 2")
}

.default_call_colours <- function() {
  c(
    "Allele 1 / Allele 1" = "#E41A1C",
    "Allele 2 / Allele 2" = "#1F4EBA",
    "Heterozygote" = "#18A535",
    "Review / no call" = "#5B5B5B",
    "NTC" = "#000000"
  )
}

.resolve_allele_2_dye <- function(kasp, allele_2_dye = NULL) {
  if (!is.null(allele_2_dye)) return(toupper(trimws(allele_2_dye)))
  dyes <- toupper(as.character(kasp$dyes))
  if ("HEX" %in% dyes) return("HEX")
  if ("VIC" %in% dyes) return("VIC")
  candidates <- setdiff(dyes, c("FAM", "ROX"))
  if (length(candidates) == 1) return(candidates)
  stop(
    "Could not infer allele_2_dye. Specify it explicitly (for example 'HEX' or 'VIC').",
    call. = FALSE
  )
}

#' Default validation settings
#'
#' Returns the default settings used for control-guided empirical clustering
#' and confidence scoring.
#'
#' @return A named list.
#' @export
kasp_validation_defaults <- function() {
  list(
    ntc_exclusion_fraction = 0.35,
    min_call_signal_fraction = 0.50,
    seed_signal_fraction = 0.65,
    seed_boundary_margin_fraction = 0.12,
    full_confidence_signal_fraction = 0.90,
    minimum_signal_confidence_cap = 40,
    min_cluster_separation = 0.10,
    local_neighbours = 3L,
    confidence_floor = 50,
    confidence_weight_boundary = 0.45,
    confidence_weight_separation = 0.35,
    confidence_weight_compactness = 0.20
  )
}
