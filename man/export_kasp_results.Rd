\name{export_kasp_results}
\alias{export_kasp_results}
\title{Export KASP analysis results}
\usage{
export_kasp_results(result, xlsx_output, html_output = NULL,
  include_coordinates = FALSE, overwrite = TRUE)
}
\arguments{
\item{result}{A \code{kasp_result} object.}
\item{xlsx_output}{Output XLSX path.}
\item{html_output}{Optional interactive HTML path.}
\item{include_coordinates}{Include the technical coordinates sheet.}
\item{overwrite}{Overwrite an existing XLSX file.}
}
\value{Invisibly returns output paths.}
\description{Exports genotype calls, controls, cluster model, settings and optional coordinates to Excel and optionally saves the interactive graph.}
