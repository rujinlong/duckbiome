% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_vir.R
\name{parse_virsorter2}
\alias{parse_virsorter2}
\title{Parse VirSorter2 Output}
\usage{
parse_virsorter2(fin_virsorter2_score, fin_virsorter2_category)
}
\arguments{
\item{fin_virsorter2_score}{A character string specifying the path to
the VirSorter2 score file.}

\item{fin_virsorter2_category}{A character string specifying the path to
the VirSorter2 category file.}
}
\value{
A data frame containing the contig ID, VirSorter2 score, group,
  and category for each contig.
}
\description{
This function parses VirSorter2 output files, including the score and
category files, to extract relevant information about viral contigs.
The function combines information from both files into a single data frame
for comprehensive analysis.
}
\details{
The function performs the following steps:
\enumerate{
  \item Reads the VirSorter2 score file and renames columns for consistency.
  \item Reads the VirSorter2 category file and renames columns.
  \item Merges the score and category information on the contig ID.
  \item Sorts the final data frame by contig ID for easier reference.
}
}
\examples{
\dontrun{
# Example usage:
# Assuming 'virsorter2_score.txt' and 'virsorter2_category.txt'
# are the VirSorter2 output files
df_virsorter2 <- parse_virsorter2(
  "virsorter2_score.txt", "virsorter2_category.txt"
)
}
}
