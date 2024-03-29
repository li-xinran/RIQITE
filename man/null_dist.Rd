% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RI_Quantile_20220919.R
\name{null_dist}
\alias{null_dist}
\title{Generate randomization dist of the rank score stat}
\usage{
null_dist(n, m, method.list = NULL, score = NULL, nperm = 10^5, Z.perm = NULL)
}
\arguments{
\item{score}{The list of possible scores of the units (under a null,
invariant under treatment assignment).  If null will generate from the
passed method.list and given n, m.}

\item{Z.perm}{Matrix of possible treatment assignments.  NULL will generate
complete rand with nperm combinations.}
}
\description{
Generate the null distribution of the given test statistic for an
experiment with m treated out of n units.
}
