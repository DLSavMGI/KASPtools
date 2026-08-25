\name{kasp_cartesian}
\alias{kasp_cartesian}
\alias{analyse_kasp}
\title{Allelic discrimination analysis}
\usage{
kasp_cartesian(kasp, plateau_start, plateau_end,
  allele_1_control_wells = NULL, allele_2_control_wells = NULL,
  heterozygote_control_wells = NULL, ntc_wells = NULL,
  allele_1_name = "Allele1", allele_2_name = "Allele2",
  heterozygote_order = c("allele_1", "allele_2"), genotype_separator = "",
  normalization = c("delta", "relative_delta", "rox"), allele_1_dye = "FAM",
  allele_2_dye = NULL, rox_dye = "ROX", min_plateau_points = 8L,
  validation_settings = list(), call_colours = NULL, html_output = NULL)

analyse_kasp(...)
}
\arguments{
\item{kasp}{A \code{kasp_run} object.}
\item{plateau_start, plateau_end}{Acquisition limits for plateau baseline calculation.}
\item{allele_1_control_wells, allele_2_control_wells, heterozygote_control_wells, ntc_wells}{Optional control-well positions.}
\item{allele_1_name, allele_2_name}{Allele labels used in exported genotype strings.}
\item{heterozygote_order}{Order of allele labels in heterozygotes.}
\item{genotype_separator}{Separator between allele labels.}
\item{normalization}{One of \code{delta}, \code{relative_delta}, or \code{rox}.}
\item{allele_1_dye, allele_2_dye, rox_dye}{Fluorescence channel names.}
\item{min_plateau_points}{Minimum finite plateau observations per well and allele channel.}
\item{validation_settings}{Named list overriding default validation settings.}
\item{call_colours}{Optional named call-colour vector.}
\item{html_output}{Optional HTML output path.}
\item{...}{Arguments passed to \code{kasp_cartesian}.}
}
\value{A \code{kasp_result} object.}
\description{Performs plateau correction, normalization, control-guided coordinate transformation, empirical clustering, confidence scoring and interactive plotting.}
