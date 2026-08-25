#' Plot fluorescence kinetics
#'
#' Plots every acquisition point retained in a `kasp_run`, including terminal
#' low-temperature acquisitions from DT96 RDML files.
#'
#' @param kasp A `kasp_run` object.
#' @param dyes Optional vector of dye names. By default all detected dyes are shown.
#' @param show_median Draw the median trajectory across wells for each dye.
#' @param show_block_boundaries Draw vertical lines between acquisition blocks when available.
#'
#' @return A `ggplot2` object.
#' @export
plot_kasp_kinetics <- function(
  kasp,
  dyes = NULL,
  show_median = TRUE,
  show_block_boundaries = TRUE
) {
  .validate_kasp_run(kasp)
  dat <- kasp$data
  available_dyes <- unique(as.character(dat$dye))

  if (is.null(dyes)) {
    preferred_order <- c("FAM", "HEX", "VIC", "ROX")
    dyes <- c(
      intersect(preferred_order, available_dyes),
      setdiff(available_dyes, preferred_order)
    )
  } else {
    dyes <- toupper(trimws(as.character(dyes)))
    missing_dyes <- setdiff(dyes, available_dyes)
    if (length(missing_dyes) > 0) {
      warning("Dyes not found: ", paste(missing_dyes, collapse = ", "), call. = FALSE)
    }
    dyes <- intersect(dyes, available_dyes)
  }
  if (length(dyes) == 0) stop("No dyes remain for plotting.", call. = FALSE)

  plot_data <- dat |>
    dplyr::filter(dye %in% dyes) |>
    dplyr::mutate(dye = factor(dye, levels = dyes))

  block_boundaries <- tibble::tibble(x = numeric())
  if (all(c("acquisition", "cycle_block") %in% names(plot_data))) {
    block_data <- plot_data |>
      dplyr::distinct(acquisition, cycle_block) |>
      dplyr::arrange(acquisition)
    if (nrow(block_data) > 1) {
      block_boundaries <- block_data |>
        dplyr::mutate(previous_block = dplyr::lag(cycle_block)) |>
        dplyr::filter(!is.na(previous_block), cycle_block != previous_block) |>
        dplyr::transmute(x = acquisition - 0.5)
    }
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = acquisition, y = fluorescence, group = well)
  ) +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.45)

  if (isTRUE(show_median)) {
    p <- p + ggplot2::stat_summary(
      ggplot2::aes(group = 1),
      fun = stats::median,
      geom = "line",
      linewidth = 1.15,
      na.rm = TRUE
    )
  }

  if (isTRUE(show_block_boundaries) && nrow(block_boundaries) > 0) {
    p <- p + ggplot2::geom_vline(
      data = block_boundaries,
      ggplot2::aes(xintercept = x),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.4
    )
  }

  p +
    ggplot2::facet_wrap(~dye, scales = "free_y", ncol = 1) +
    ggplot2::scale_x_continuous(
      breaks = seq(
        min(plot_data$acquisition, na.rm = TRUE),
        max(plot_data$acquisition, na.rm = TRUE),
        by = 5
      ),
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(
      title = paste0("Fluorescence kinetics — ", kasp$protocol_name),
      x = "Acquisition",
      y = "Fluorescence"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}
