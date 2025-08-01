% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_vir.R
\name{parse_vcontact3}
\alias{parse_vcontact3}
\title{Parse vContact3 Output}
\usage{
parse_vcontact3(fin_vcontact3)
}
\arguments{
\item{fin_vcontact3}{Character string specifying the path to the
vContact3 `final_assignments.csv` file.}
}
\value{
A data frame containing parsed vContact3 information with
  VC IDs, classification status, reference presence, and novel
  genus flags.
}
\description{
This function parses the vContact3 output file, specifically the
`final_assignments.csv` file. It processes genus predictions,
classifies contigs, and identifies novel genera with or without
reference genomes.
}
\examples{
\dontrun{
# Assuming 'final_assignments.csv' is the vContact3 output file
df_vcontact3 <- parse_vcontact3("final_assignments.csv")
}
}
