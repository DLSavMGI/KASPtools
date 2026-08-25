\name{import_kasp_xml}
\alias{import_kasp_xml}
\title{Import DT96 RDML/XML fluorescence data}
\usage{
import_kasp_xml(file, protocol_to_keep = 1L, protocol_name, columns = NULL,
  output_rds = NULL)
}
\arguments{
\item{file}{Path to an RDML/XML file.}
\item{protocol_to_keep}{Automatically detected assay-block index.}
\item{protocol_name}{User-facing protocol name stored in the object.}
\item{columns}{Optional plate columns to select explicitly.}
\item{output_rds}{Optional RDS output path.}
}
\value{A \code{kasp_run} object.}
\description{Imports all fluorescence channels and acquisition points from a DT96 RDML/XML export.}
