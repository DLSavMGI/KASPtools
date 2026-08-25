#' Allelic discrimination analysis for KASP/KASP-like data
#'
#' Applies per-well plateau correction, constructs endpoint coordinates,
#' optionally normalizes by ROX or the per-channel plateau, transforms the
#' coordinates into a control-guided basis, learns empirical genotype clusters,
#' calculates confidence scores, and creates an interactive Plotly scatter plot.
#'
#' @param kasp A `kasp_run` object.
#' @param plateau_start First acquisition used for the plateau baseline.
#' @param plateau_end Last acquisition used for the plateau baseline.
#' @param allele_1_control_wells Optional wells for allele-1 homozygous controls.
#' @param allele_2_control_wells Optional wells for allele-2 homozygous controls.
#' @param heterozygote_control_wells Optional wells for heterozygous controls.
#' @param ntc_wells Optional NTC wells. QuantStudio Task annotations are used when available.
#' @param allele_1_name Genotype label for allele 1.
#' @param allele_2_name Genotype label for allele 2.
#' @param heterozygote_order Order of allele labels in the heterozygote genotype string.
#' @param genotype_separator Separator used between allele labels.
#' @param normalization One of `"delta"`, `"relative_delta"`, or `"rox"`.
#' @param allele_1_dye Dye representing allele 1.
#' @param allele_2_dye Dye representing allele 2. If NULL, HEX is preferred, then VIC.
#' @param rox_dye Passive-reference dye used when `normalization = "rox"`.
#' @param min_plateau_points Minimum number of finite values required in the plateau window.
#' @param validation_settings Named list overriding values from `kasp_validation_defaults()`.
#' @param call_colours Optional named vector overriding default plot colours.
#' @param html_output Optional path for saving the interactive plot as HTML.
#'
#' @return A `kasp_result` object containing the plot, calls, coordinates and model diagnostics.
#' @export
kasp_cartesian <- function(
  kasp,
  plateau_start,
  plateau_end,
  allele_1_control_wells = NULL,
  allele_2_control_wells = NULL,
  heterozygote_control_wells = NULL,
  ntc_wells = NULL,
  allele_1_name = "Allele1",
  allele_2_name = "Allele2",
  heterozygote_order = c("allele_1", "allele_2"),
  genotype_separator = "",
  normalization = c("delta", "relative_delta", "rox"),
  allele_1_dye = "FAM",
  allele_2_dye = NULL,
  rox_dye = "ROX",
  min_plateau_points = 8L,
  validation_settings = list(),
  call_colours = NULL,
  html_output = NULL
) {
  .validate_kasp_run(kasp)
  normalization <- match.arg(normalization)
  allele_1_dye <- toupper(trimws(as.character(allele_1_dye)))
  allele_2_dye <- .resolve_allele_2_dye(kasp, allele_2_dye)
  rox_dye <- toupper(trimws(as.character(rox_dye)))
  if (allele_1_dye == allele_2_dye) stop("allele_1_dye and allele_2_dye must differ.", call. = FALSE)

  settings <- utils::modifyList(kasp_validation_defaults(), validation_settings)
  numeric_settings <- c(
    "ntc_exclusion_fraction", "min_call_signal_fraction", "seed_signal_fraction",
    "seed_boundary_margin_fraction", "full_confidence_signal_fraction",
    "minimum_signal_confidence_cap", "min_cluster_separation", "local_neighbours",
    "confidence_floor", "confidence_weight_boundary", "confidence_weight_separation",
    "confidence_weight_compactness"
  )
  missing_settings <- setdiff(numeric_settings, names(settings))
  if (length(missing_settings) > 0) {
    stop("Missing validation settings: ", paste(missing_settings, collapse = ", "), call. = FALSE)
  }
  fraction_names <- c(
    "ntc_exclusion_fraction", "min_call_signal_fraction", "seed_signal_fraction",
    "seed_boundary_margin_fraction", "full_confidence_signal_fraction", "min_cluster_separation"
  )
  fraction_values <- unlist(settings[fraction_names], use.names = TRUE)
  if (any(!is.finite(fraction_values)) || any(fraction_values < 0)) {
    stop("Validation fractions must be finite and non-negative.", call. = FALSE)
  }
  if (settings$min_cluster_separation >= 1) {
    stop("min_cluster_separation must be smaller than 1.", call. = FALSE)
  }
  if (settings$full_confidence_signal_fraction <= settings$min_call_signal_fraction) {
    stop("full_confidence_signal_fraction must exceed min_call_signal_fraction.", call. = FALSE)
  }
  if (!is.finite(settings$minimum_signal_confidence_cap) ||
      settings$minimum_signal_confidence_cap < 0 || settings$minimum_signal_confidence_cap > 100) {
    stop("minimum_signal_confidence_cap must be between 0 and 100.", call. = FALSE)
  }
  if (!is.finite(settings$confidence_floor) || settings$confidence_floor < 0 || settings$confidence_floor > 100) {
    stop("confidence_floor must be between 0 and 100.", call. = FALSE)
  }
  if (!is.finite(settings$local_neighbours) || settings$local_neighbours < 1) {
    stop("local_neighbours must be at least 1.", call. = FALSE)
  }
  weights <- c(
    settings$confidence_weight_boundary,
    settings$confidence_weight_separation,
    settings$confidence_weight_compactness
  )
  if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("Confidence weights must be finite, non-negative, and have a positive sum.", call. = FALSE)
  }
  weights <- weights / sum(weights)

  plate_map <- kasp$plate_map
  if (!"task" %in% names(plate_map)) plate_map$task <- "UNKNOWN"
  if (!"omit" %in% names(plate_map)) plate_map$omit <- FALSE
  plate_map$well_position <- toupper(trimws(as.character(plate_map$well_position)))
  plate_map$task <- toupper(trimws(as.character(plate_map$task)))
  plate_map$task[is.na(plate_map$task) | plate_map$task == ""] <- "UNKNOWN"
  plate_map$omit[is.na(plate_map$omit)] <- FALSE
  plate_map <- plate_map[!plate_map$omit, , drop = FALSE]

  a1_wells <- .normalize_wells(allele_1_control_wells)
  a2_wells <- .normalize_wells(allele_2_control_wells)
  het_wells <- .normalize_wells(heterozygote_control_wells)
  n_wells <- .normalize_wells(ntc_wells)
  supplied_control_wells <- c(a1_wells, a2_wells, het_wells, n_wells)
  duplicate_controls <- unique(supplied_control_wells[duplicated(supplied_control_wells)])
  if (length(duplicate_controls) > 0) {
    stop("Control wells were assigned to multiple classes: ", paste(duplicate_controls, collapse = ", "), call. = FALSE)
  }
  missing_control_wells <- setdiff(supplied_control_wells, plate_map$well_position)
  if (length(missing_control_wells) > 0) {
    stop("Control wells are absent from the selected assay: ", paste(missing_control_wells, collapse = ", "), call. = FALSE)
  }

  if (length(a1_wells) > 0) plate_map$task[plate_map$well_position %in% a1_wells] <- "PC_ALLELE_1"
  if (length(a2_wells) > 0) plate_map$task[plate_map$well_position %in% a2_wells] <- "PC_ALLELE_2"
  if (length(het_wells) > 0) plate_map$task[plate_map$well_position %in% het_wells] <- "PC_ALLELE_BOTH"
  if (length(n_wells) > 0) plate_map$task[plate_map$well_position %in% n_wells] <- "NTC"

  available_dyes <- toupper(as.character(kasp$dyes))
  required_dyes <- c(allele_1_dye, allele_2_dye)
  if (normalization == "rox") required_dyes <- c(required_dyes, rox_dye)
  missing_dyes <- setdiff(required_dyes, available_dyes)
  if (length(missing_dyes) > 0) {
    stop(
      "Normalization '", normalization, "' requires missing dye(s): ",
      paste(missing_dyes, collapse = ", "),
      call. = FALSE
    )
  }

  dat <- kasp$data
  dat$dye <- toupper(trimws(as.character(dat$dye)))
  dat <- dat[dat$well %in% plate_map$well, , drop = FALSE]

  if (!is.finite(plateau_start) || !is.finite(plateau_end) || plateau_start > plateau_end) {
    stop("plateau_start and plateau_end must define a valid increasing interval.", call. = FALSE)
  }
  plateau_requested <- seq.int(as.integer(plateau_start), as.integer(plateau_end))
  available_acquisitions <- sort(unique(dat$acquisition[is.finite(dat$acquisition)]))
  plateau_points <- intersect(plateau_requested, available_acquisitions)
  if (length(plateau_points) < as.integer(min_plateau_points)) {
    stop(
      "Plateau window contains only ", length(plateau_points),
      " available acquisitions; minimum is ", min_plateau_points, ".",
      call. = FALSE
    )
  }

  baseline_long <- dat |>
    dplyr::filter(
      dye %in% c(allele_1_dye, allele_2_dye),
      acquisition %in% plateau_points
    ) |>
    dplyr::group_by(well, dye) |>
    dplyr::summarise(
      n_plateau = sum(is.finite(fluorescence)),
      plateau = if (sum(is.finite(fluorescence)) >= as.integer(min_plateau_points)) {
        stats::median(fluorescence, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  baseline_wide <- baseline_long |>
    dplyr::select(well, dye, plateau) |>
    tidyr::pivot_wider(
      names_from = dye,
      values_from = plateau,
      names_glue = "{dye}_plateau"
    )

  endpoint_acquisition <- max(dat$acquisition, na.rm = TRUE)
  endpoint_long <- dat |>
    dplyr::filter(acquisition == endpoint_acquisition, dye %in% required_dyes) |>
    dplyr::select(well, dye, fluorescence)
  endpoint_wide <- endpoint_long |>
    tidyr::pivot_wider(
      names_from = dye,
      values_from = fluorescence,
      names_glue = "{dye}_endpoint"
    )

  plot_base <- plate_map |>
    dplyr::left_join(baseline_wide, by = "well") |>
    dplyr::left_join(endpoint_wide, by = "well")

  a1_plateau_column <- paste0(allele_1_dye, "_plateau")
  a2_plateau_column <- paste0(allele_2_dye, "_plateau")
  a1_endpoint_column <- paste0(allele_1_dye, "_endpoint")
  a2_endpoint_column <- paste0(allele_2_dye, "_endpoint")
  required_coordinate_columns <- c(
    a1_plateau_column, a2_plateau_column,
    a1_endpoint_column, a2_endpoint_column
  )
  missing_coordinate_columns <- setdiff(required_coordinate_columns, names(plot_base))
  if (length(missing_coordinate_columns) > 0) {
    stop("Missing fluorescence values needed for coordinates: ", paste(missing_coordinate_columns, collapse = ", "), call. = FALSE)
  }

  plot_base$A1_Plateau <- plot_base[[a1_plateau_column]]
  plot_base$A2_Plateau <- plot_base[[a2_plateau_column]]
  plot_base$A1_Endpoint <- plot_base[[a1_endpoint_column]]
  plot_base$A2_Endpoint <- plot_base[[a2_endpoint_column]]
  plot_base$A1_Delta <- plot_base$A1_Endpoint - plot_base$A1_Plateau
  plot_base$A2_Delta <- plot_base$A2_Endpoint - plot_base$A2_Plateau

  if (normalization == "delta") {
    plot_base$X <- plot_base$A1_Delta
    plot_base$Y <- plot_base$A2_Delta
  } else if (normalization == "relative_delta") {
    plot_base$X <- plot_base$A1_Delta / plot_base$A1_Plateau
    plot_base$Y <- plot_base$A2_Delta / plot_base$A2_Plateau
  } else {
    rox_endpoint_column <- paste0(rox_dye, "_endpoint")
    if (!rox_endpoint_column %in% names(plot_base)) {
      stop("ROX endpoint data are missing.", call. = FALSE)
    }
    plot_base$ROX_Endpoint <- plot_base[[rox_endpoint_column]]
    plot_base$X <- plot_base$A1_Delta / plot_base$ROX_Endpoint
    plot_base$Y <- plot_base$A2_Delta / plot_base$ROX_Endpoint
  }

  n_before <- nrow(plot_base)
  plot_base <- plot_base[is.finite(plot_base$X) & is.finite(plot_base$Y), , drop = FALSE]
  n_removed <- n_before - nrow(plot_base)
  if (n_removed > 0) warning(n_removed, " wells were removed because coordinates were not finite.", call. = FALSE)
  if (nrow(plot_base) == 0) stop("No finite endpoint coordinates remain.", call. = FALSE)
  plot_base$Hover <- .make_hover_text(plot_base$well_position, plot_base$sample_name)

  control_task_to_call <- .control_task_to_call()
  genotype_calls <- .genotype_calls()
  required_control_tasks <- names(control_task_to_call)
  missing_control_tasks <- setdiff(required_control_tasks, unique(plot_base$task))
  if (length(missing_control_tasks) > 0) {
    stop(
      "Required control classes are missing: ", paste(missing_control_tasks, collapse = ", "),
      ". Supply control well positions if the instrument export does not contain Task annotations.",
      call. = FALSE
    )
  }

  controls_raw <- plot_base[plot_base$task %in% required_control_tasks, , drop = FALSE]
  controls_raw$Cluster_Call <- unname(control_task_to_call[controls_raw$task])
  samples_raw <- plot_base[plot_base$task == "UNKNOWN", , drop = FALSE]
  ntc_raw <- plot_base[plot_base$task == "NTC", , drop = FALSE]
  if (nrow(samples_raw) == 0) stop("No UNKNOWN samples were found.", call. = FALSE)
  if (nrow(ntc_raw) == 0) stop("At least one NTC is required.", call. = FALSE)

  .raw_center <- function(x) {
    c(
      X = stats::median(x$X, na.rm = TRUE),
      Y = stats::median(x$Y, na.rm = TRUE)
    )
  }
  ntc_center_raw <- .raw_center(ntc_raw)
  a1_center_raw <- .raw_center(controls_raw[controls_raw$Cluster_Call == "Allele 1 / Allele 1", , drop = FALSE])
  a2_center_raw <- .raw_center(controls_raw[controls_raw$Cluster_Call == "Allele 2 / Allele 2", , drop = FALSE])

  basis_matrix <- cbind(a1_center_raw - ntc_center_raw, a2_center_raw - ntc_center_raw)
  if (!all(is.finite(basis_matrix)) || abs(det(basis_matrix)) < 1e-10) {
    stop("Allele-control directions are nearly collinear; control-guided coordinates cannot be constructed.", call. = FALSE)
  }

  components <- .transform_to_control_basis(
    x = plot_base$X,
    y = plot_base$Y,
    ntc_x = ntc_center_raw[["X"]],
    ntc_y = ntc_center_raw[["Y"]],
    basis_matrix = basis_matrix
  )
  plot_base <- dplyr::bind_cols(plot_base, components)
  plot_base$Signal_Strength <- sqrt(plot_base$A1_Component^2 + plot_base$A2_Component^2)
  plot_base$Component_Sum <- plot_base$A1_Component + plot_base$A2_Component
  plot_base$Allele1_Fraction <- ifelse(
    is.finite(plot_base$Component_Sum) & plot_base$Component_Sum > 0,
    plot_base$A1_Component / plot_base$Component_Sum,
    NA_real_
  )

  controls_df <- plot_base[plot_base$task %in% required_control_tasks, , drop = FALSE]
  controls_df$Cluster_Call <- unname(control_task_to_call[controls_df$task])
  samples_df <- plot_base[plot_base$task == "UNKNOWN", , drop = FALSE]
  ntc_df <- plot_base[plot_base$task == "NTC", , drop = FALSE]
  ntc_df$Cluster_Call <- "NTC"

  control_reference <- controls_df |>
    dplyr::group_by(Cluster_Call) |>
    dplyr::summarise(
      Control_Balance = stats::median(Allele1_Fraction, na.rm = TRUE),
      Control_Signal = stats::median(Signal_Strength, na.rm = TRUE),
      n_controls = dplyr::n(),
      .groups = "drop"
    )

  .control_value <- function(cluster_name, column) {
    z <- control_reference[control_reference$Cluster_Call == cluster_name, column, drop = TRUE]
    if (length(z) != 1L || !is.finite(z)) {
      stop("Could not determine control centre for: ", cluster_name, call. = FALSE)
    }
    z
  }
  balance_a1 <- .control_value("Allele 1 / Allele 1", "Control_Balance")
  balance_hetero <- .control_value("Heterozygote", "Control_Balance")
  balance_a2 <- .control_value("Allele 2 / Allele 2", "Control_Balance")
  if (!(balance_a1 > balance_hetero && balance_hetero > balance_a2)) {
    stop("Control centres are not ordered Allele 1 -> Heterozygote -> Allele 2. Check dyes and control wells.", call. = FALSE)
  }

  positive_control_signals <- control_reference$Control_Signal[
    is.finite(control_reference$Control_Signal) & control_reference$Control_Signal > 0
  ]
  if (length(positive_control_signals) == 0) stop("Positive-control signal reference could not be determined.", call. = FALSE)
  positive_signal_reference <- min(positive_control_signals)
  ntc_exclusion_radius <- settings$ntc_exclusion_fraction * positive_signal_reference
  min_signal_for_call <- settings$min_call_signal_fraction * positive_signal_reference
  min_signal_for_seed <- settings$seed_signal_fraction * positive_signal_reference
  full_signal_for_confidence <- settings$full_confidence_signal_fraction * positive_signal_reference

  initial_boundary_a1_het <- mean(c(balance_a1, balance_hetero))
  initial_boundary_het_a2 <- mean(c(balance_hetero, balance_a2))
  seed_margin_a1_het <- settings$seed_boundary_margin_fraction * (balance_a1 - balance_hetero)
  seed_margin_het_a2 <- settings$seed_boundary_margin_fraction * (balance_hetero - balance_a2)

  samples_seeded <- samples_df |>
    dplyr::mutate(
      Seed_Call = dplyr::case_when(
        Signal_Strength >= min_signal_for_seed &
          Allele1_Fraction >= initial_boundary_a1_het + seed_margin_a1_het ~ "Allele 1 / Allele 1",
        Signal_Strength >= min_signal_for_seed &
          Allele1_Fraction <= initial_boundary_het_a2 - seed_margin_het_a2 ~ "Allele 2 / Allele 2",
        Signal_Strength >= min_signal_for_seed &
          Allele1_Fraction < initial_boundary_a1_het - seed_margin_a1_het &
          Allele1_Fraction > initial_boundary_het_a2 + seed_margin_het_a2 ~ "Heterozygote",
        TRUE ~ NA_character_
      )
    )

  training_points <- dplyr::bind_rows(
    controls_df |>
      dplyr::transmute(Cluster_Call, Allele1_Fraction),
    samples_seeded |>
      dplyr::filter(!is.na(Seed_Call)) |>
      dplyr::transmute(Cluster_Call = Seed_Call, Allele1_Fraction)
  )

  cluster_model <- training_points |>
    dplyr::group_by(Cluster_Call) |>
    dplyr::summarise(
      Balance_Center = stats::median(Allele1_Fraction, na.rm = TRUE),
      n_training = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(match(Cluster_Call, genotype_calls))
  if (nrow(cluster_model) != 3L || any(!is.finite(cluster_model$Balance_Center))) {
    stop("Could not form all three empirical genotype clusters.", call. = FALSE)
  }

  .cluster_balance <- function(name) {
    z <- cluster_model$Balance_Center[cluster_model$Cluster_Call == name]
    if (length(z) != 1L || !is.finite(z)) stop("Missing empirical cluster: ", name, call. = FALSE)
    z
  }
  emp_balance_a1 <- .cluster_balance("Allele 1 / Allele 1")
  emp_balance_hetero <- .cluster_balance("Heterozygote")
  emp_balance_a2 <- .cluster_balance("Allele 2 / Allele 2")
  if (!(emp_balance_a1 > emp_balance_hetero && emp_balance_hetero > emp_balance_a2)) {
    stop("Empirical clusters are not ordered correctly. Check the plateau window and seed settings.", call. = FALSE)
  }
  boundary_a1_het <- mean(c(emp_balance_a1, emp_balance_hetero))
  boundary_het_a2 <- mean(c(emp_balance_hetero, emp_balance_a2))

  n_samples <- nrow(samples_df)
  distance_matrix <- matrix(NA_real_, nrow = n_samples, ncol = nrow(cluster_model))
  colnames(distance_matrix) <- cluster_model$Cluster_Call
  valid_balance_index <- which(is.finite(samples_df$Allele1_Fraction))
  if (length(valid_balance_index) > 0) {
    distance_matrix[valid_balance_index, ] <- vapply(
      cluster_model$Balance_Center,
      function(center) abs(samples_df$Allele1_Fraction[valid_balance_index] - center),
      numeric(length(valid_balance_index))
    )
  }

  nearest_cluster <- rep(NA_character_, n_samples)
  nearest_distance <- rep(NA_real_, n_samples)
  second_distance <- rep(NA_real_, n_samples)
  relative_separation <- rep(NA_real_, n_samples)
  if (length(valid_balance_index) > 0) {
    valid_distances <- distance_matrix[valid_balance_index, , drop = FALSE]
    nearest_index <- apply(valid_distances, 1, which.min)
    nearest_cluster[valid_balance_index] <- cluster_model$Cluster_Call[nearest_index]
    nearest_distance[valid_balance_index] <- valid_distances[cbind(seq_along(valid_balance_index), nearest_index)]
    second_distance[valid_balance_index] <- apply(valid_distances, 1, function(x) sort(x)[2])
    relative_separation[valid_balance_index] <- 1 - nearest_distance[valid_balance_index] /
      pmax(second_distance[valid_balance_index], .Machine$double.eps)
  }

  samples_validated <- samples_df |>
    dplyr::mutate(
      Empirical_Cluster = nearest_cluster,
      Distance_To_Nearest_Cluster = nearest_distance,
      Distance_To_Second_Cluster = second_distance,
      Relative_Cluster_Separation = relative_separation,
      Flag_Invalid_Geometry = !is.finite(Allele1_Fraction) |
        !is.finite(Signal_Strength) | is.na(Empirical_Cluster),
      Flag_Near_NTC = is.finite(Signal_Strength) & Signal_Strength < ntc_exclusion_radius,
      Flag_Low_Signal = is.finite(Signal_Strength) & Signal_Strength < min_signal_for_call,
      Flag_Ambiguous_Between_Clusters = is.finite(Relative_Cluster_Separation) &
        Relative_Cluster_Separation < settings$min_cluster_separation,
      Is_Valid_Call = !Flag_Invalid_Geometry & !Flag_Near_NTC &
        !Flag_Low_Signal & !Flag_Ambiguous_Between_Clusters,
      Final_Cluster_Call = dplyr::if_else(Is_Valid_Call, Empirical_Cluster, "Review / no call")
    )

  samples_validated <- .calculate_local_compactness(samples_validated, k = settings$local_neighbours)

  samples_validated <- samples_validated |>
    dplyr::mutate(
      Score_Boundary = dplyr::case_when(
        Empirical_Cluster == "Allele 1 / Allele 1" ~ .clamp01(
          (Allele1_Fraction - boundary_a1_het) /
            pmax(emp_balance_a1 - boundary_a1_het, .Machine$double.eps)
        ),
        Empirical_Cluster == "Allele 2 / Allele 2" ~ .clamp01(
          (boundary_het_a2 - Allele1_Fraction) /
            pmax(boundary_het_a2 - emp_balance_a2, .Machine$double.eps)
        ),
        Empirical_Cluster == "Heterozygote" ~ .clamp01(
          pmin(
            (Allele1_Fraction - boundary_het_a2) /
              pmax(emp_balance_hetero - boundary_het_a2, .Machine$double.eps),
            (boundary_a1_het - Allele1_Fraction) /
              pmax(boundary_a1_het - emp_balance_hetero, .Machine$double.eps)
          )
        ),
        TRUE ~ 0
      ),
      Score_Separation = .clamp01(
        (Relative_Cluster_Separation - settings$min_cluster_separation) /
          pmax(1 - settings$min_cluster_separation, .Machine$double.eps)
      ),
      Cluster_Quality_Score = .clamp01(
        weights[1] * Score_Boundary +
          weights[2] * Score_Separation +
          weights[3] * Score_Local_Compactness
      ),
      Cluster_Confidence_Score = settings$confidence_floor +
        (100 - settings$confidence_floor) * Cluster_Quality_Score,
      Score_Signal_Clearance = .clamp01(
        (Signal_Strength - min_signal_for_call) /
          pmax(full_signal_for_confidence - min_signal_for_call, .Machine$double.eps)
      ),
      Signal_Confidence_Cap = settings$minimum_signal_confidence_cap +
        (100 - settings$minimum_signal_confidence_cap) * Score_Signal_Clearance,
      Call_Confidence_Score = dplyr::case_when(
        Is_Valid_Call ~ round(pmin(Cluster_Confidence_Score, Signal_Confidence_Cap), 1),
        TRUE ~ 0
      ),
      Confidence_Category = dplyr::case_when(
        !Is_Valid_Call ~ "Review",
        Call_Confidence_Score >= 80 ~ "High",
        Call_Confidence_Score >= 65 ~ "Moderate",
        TRUE ~ "Low"
      ),
      Review_Note = dplyr::case_when(
        Flag_Invalid_Geometry ~ "Invalid signal geometry",
        Flag_Near_NTC ~ "Near NTC",
        Flag_Low_Signal ~ "Low signal",
        Flag_Ambiguous_Between_Clusters ~ "Between genotype clusters",
        TRUE ~ ""
      )
    )

  allele_labels <- c(allele_1 = allele_1_name, allele_2 = allele_2_name)
  if (length(heterozygote_order) != 2L || !setequal(heterozygote_order, names(allele_labels))) {
    stop("heterozygote_order must contain allele_1 and allele_2 exactly once.", call. = FALSE)
  }
  genotype_a1 <- paste0(allele_1_name, genotype_separator, allele_1_name)
  genotype_a2 <- paste0(allele_2_name, genotype_separator, allele_2_name)
  genotype_hetero <- paste0(
    allele_labels[heterozygote_order[1]], genotype_separator,
    allele_labels[heterozygote_order[2]]
  )
  call_to_genotype <- c(
    "Allele 1 / Allele 1" = genotype_a1,
    "Allele 2 / Allele 2" = genotype_a2,
    "Heterozygote" = genotype_hetero,
    "Review / no call" = "Undetermined",
    "NTC" = "NTC"
  )

  samples_validated$Genotype <- unname(call_to_genotype[samples_validated$Final_Cluster_Call])
  samples_validated$Plotly_Hover <- ifelse(
    samples_validated$Is_Valid_Call,
    paste0(
      samples_validated$Hover,
      "<br><b>Genotype:</b> ", samples_validated$Genotype,
      "<br><b>Confidence:</b> ", samples_validated$Call_Confidence_Score, "%"
    ),
    paste0(
      samples_validated$Hover,
      "<br><b>Genotype:</b> Undetermined",
      "<br><b>Review:</b> ", samples_validated$Review_Note
    )
  )

  controls_processed <- controls_df
  controls_processed$Final_Cluster_Call <- controls_processed$Cluster_Call
  controls_processed$Call_Confidence_Score <- NA_real_
  controls_processed$Confidence_Category <- "Control"
  controls_processed$Review_Note <- ""
  controls_processed$Is_Valid_Call <- TRUE
  controls_processed$Genotype <- unname(call_to_genotype[controls_processed$Cluster_Call])
  controls_processed$Plotly_Hover <- paste0(
    controls_processed$Hover,
    "<br><b>Genotype:</b> ", controls_processed$Genotype,
    "<br><b>Type:</b> Positive control"
  )

  ntc_processed <- ntc_df
  ntc_processed$Final_Cluster_Call <- "NTC"
  ntc_processed$Call_Confidence_Score <- NA_real_
  ntc_processed$Confidence_Category <- "NTC"
  ntc_processed$Review_Note <- ""
  ntc_processed$Is_Valid_Call <- TRUE
  ntc_processed$Genotype <- "NTC"
  ntc_processed$Plotly_Hover <- paste0(ntc_processed$Hover, "<br><b>Type:</b> NTC")

  colours <- .default_call_colours()
  if (!is.null(call_colours)) {
    if (is.null(names(call_colours))) stop("call_colours must be a named vector.", call. = FALSE)
    colours[names(call_colours)] <- call_colours
  }

  valid_samples <- samples_validated[samples_validated$Is_Valid_Call, , drop = FALSE]
  review_samples <- samples_validated[!samples_validated$Is_Valid_Call, , drop = FALSE]
  x_axis_range <- .make_axis_range(plot_base$Y)
  y_axis_range <- .make_axis_range(plot_base$X)

  x_axis_title <- switch(
    normalization,
    delta = paste0(allele_2_dye, ": endpoint - plateau"),
    relative_delta = paste0(allele_2_dye, ": (endpoint - plateau) / plateau"),
    rox = paste0(allele_2_dye, ": (endpoint - plateau) / ", rox_dye, " endpoint")
  )
  y_axis_title <- switch(
    normalization,
    delta = paste0(allele_1_dye, ": endpoint - plateau"),
    relative_delta = paste0(allele_1_dye, ": (endpoint - plateau) / plateau"),
    rox = paste0(allele_1_dye, ": (endpoint - plateau) / ", rox_dye, " endpoint")
  )

  p <- plotly::plot_ly()
  for (call_name in genotype_calls) {
    current_samples <- valid_samples[valid_samples$Final_Cluster_Call == call_name, , drop = FALSE]
    if (nrow(current_samples) > 0) {
      p <- plotly::add_markers(
        p,
        data = current_samples,
        x = ~Y,
        y = ~X,
        text = ~Plotly_Hover,
        hoverinfo = "text",
        name = unname(call_to_genotype[call_name]),
        legendgroup = call_name,
        marker = list(size = 9, color = unname(colours[call_name]))
      )
    }
  }
  if (nrow(review_samples) > 0) {
    p <- plotly::add_markers(
      p,
      data = review_samples,
      x = ~Y,
      y = ~X,
      text = ~Plotly_Hover,
      hoverinfo = "text",
      name = "Undetermined",
      legendgroup = "Review / no call",
      marker = list(symbol = "x", size = 12, color = unname(colours["Review / no call"]))
    )
  }
  for (call_name in genotype_calls) {
    current_controls <- controls_processed[controls_processed$Final_Cluster_Call == call_name, , drop = FALSE]
    if (nrow(current_controls) > 0) {
      class_exists <- any(valid_samples$Final_Cluster_Call == call_name)
      p <- plotly::add_markers(
        p,
        data = current_controls,
        x = ~Y,
        y = ~X,
        text = ~Plotly_Hover,
        hoverinfo = "text",
        name = unname(call_to_genotype[call_name]),
        legendgroup = call_name,
        showlegend = !class_exists,
        marker = list(
          size = 16,
          color = unname(colours[call_name]),
          line = list(color = "#222222", width = 1)
        )
      )
    }
  }
  if (nrow(ntc_processed) > 0) {
    p <- plotly::add_markers(
      p,
      data = ntc_processed,
      x = ~Y,
      y = ~X,
      text = ~Plotly_Hover,
      hoverinfo = "text",
      name = "NTC",
      legendgroup = "NTC",
      marker = list(size = 14, color = unname(colours["NTC"]))
    )
  }

  p <- plotly::layout(
    p,
    title = list(
      text = paste0(
        "<b>", kasp$protocol_name, ": Allelic Discrimination Plot</b>",
        "<br><sup>Normalization: ", normalization,
        "; endpoint acquisition ", endpoint_acquisition,
        "; plateau ", min(plateau_points), "-", max(plateau_points), "</sup>"
      ),
      x = 0.5,
      xanchor = "center"
    ),
    xaxis = list(
      title = x_axis_title,
      range = x_axis_range,
      showline = TRUE,
      linecolor = "#777777",
      gridcolor = "#E5E5E5",
      zeroline = TRUE,
      zerolinecolor = "#D0D0D0"
    ),
    yaxis = list(
      title = y_axis_title,
      range = y_axis_range,
      showline = TRUE,
      linecolor = "#777777",
      gridcolor = "#E5E5E5",
      zeroline = TRUE,
      zerolinecolor = "#D0D0D0"
    ),
    legend = list(title = list(text = "Genotype"), x = 1.02, y = 1),
    hovermode = "closest",
    plot_bgcolor = "white",
    paper_bgcolor = "white",
    margin = list(l = 100, r = 190, t = 100, b = 90)
  )
  p <- plotly::config(p, displaylogo = FALSE, responsive = TRUE, scrollZoom = TRUE)

  if (!is.null(html_output)) {
    tryCatch(
      htmlwidgets::saveWidget(p, file = html_output, selfcontained = TRUE),
      error = function(e) {
        warning("Could not save self-contained HTML; saving with dependencies directory.", call. = FALSE)
        htmlwidgets::saveWidget(p, file = html_output, selfcontained = FALSE)
      }
    )
  }

  result <- list(
    protocol_name = kasp$protocol_name,
    source_file = kasp$source_file,
    source_type = kasp$source_type,
    normalization = normalization,
    allele_1_dye = allele_1_dye,
    allele_2_dye = allele_2_dye,
    rox_dye = rox_dye,
    allele_1_name = allele_1_name,
    allele_2_name = allele_2_name,
    heterozygote_order = heterozygote_order,
    genotype_separator = genotype_separator,
    plateau = plateau_points,
    endpoint_acquisition = endpoint_acquisition,
    positive_signal_reference = positive_signal_reference,
    plot = p,
    samples = samples_validated,
    controls = controls_processed,
    ntc = ntc_processed,
    coordinates = plot_base,
    cluster_model = cluster_model,
    control_reference = control_reference,
    call_to_genotype = call_to_genotype,
    settings = settings
  )
  class(result) <- c("kasp_result", "dt_kasp_result", "list")
  result
}

#' @rdname kasp_cartesian
#' @export
analyse_kasp <- function(...) kasp_cartesian(...)
