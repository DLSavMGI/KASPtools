\name{import_kasp_quantstudio}
\alias{import_kasp_quantstudio}
\title{Import QuantStudio KASP fluorescence data}
\usage{
import_kasp_quantstudio(file, assay_to_keep, protocol_name = assay_to_keep,
  dye_columns = c("FAM", "VIC", "HEX", "ROX"), include_omitted = FALSE,
  output_rds = NULL)
}
\arguments{
\item{file}{Path to a QuantStudio XLSX export.}
\item{assay_to_keep}{Exact SNP assay name.}
\item{protocol_name}{User-facing protocol name.}
\item{dye_columns}{Candidate fluorescence channels to import.}
\item{include_omitted}{Retain wells marked Omit.}
\item{output_rds}{Optional RDS output path.}
}
\value{A \code{kasp_run} object.}
\description{Reads Results and Multicomponent Data and converts one assay to the common KASPtools data structure.}
