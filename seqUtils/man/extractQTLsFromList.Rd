% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qtl_postProcess.R
\name{extractQTLsFromList}
\alias{extractQTLsFromList}
\title{Extract all QTLs at a specific FDR level from a list of min pvalues by condition}
\usage{
extractQTLsFromList(min_pvalue_list, fdr_cutoff = 0.1)
}
\arguments{
\item{min_pvalue_list}{List of QTLs per condition}

\item{fdr_cutoff}{}
}
\value{
Data frame of QTLs
}
\description{
Multiple variants per gene are sorted by p-value
}
