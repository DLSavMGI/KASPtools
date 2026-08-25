\name{plot_kasp_kinetics}
\alias{plot_kasp_kinetics}
\title{Plot fluorescence kinetics}
\usage{
plot_kasp_kinetics(kasp, dyes = NULL, show_median = TRUE,
  show_block_boundaries = TRUE)
}
\arguments{
\item{kasp}{A \code{kasp_run} object.}
\item{dyes}{Optional vector of dye names.}
\item{show_median}{Show the median trajectory.}
\item{show_block_boundaries}{Show boundaries between acquisition blocks when known.}
}
\value{A ggplot2 object.}
\description{Plots all retained fluorescence acquisition points for manual plateau-window selection.}
