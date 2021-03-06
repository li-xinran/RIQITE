% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perminf_code.R
\name{find_ci}
\alias{find_ci}
\alias{find_ci_stephenson}
\title{Confidence intervals for maximal (minimal) effects}
\usage{
find_ci(
  data,
  coin_test,
  formula,
  tr_level = NULL,
  alternative = c("greater", "less"),
  distribution = "asymptotic",
  ci_level = 0.9,
  step_fraction = 0.001,
  precision = 0.001,
  verbosity = 0,
  round_digits = 3,
  init_tau0 = NULL,
  tau_grid = NULL,
  pvalue_ci = TRUE,
  ...
)

find_ci_stephenson(
  data,
  formula,
  tr_level = NULL,
  distribution = "asymptotic",
  alternative = "less",
  tail = c("right", "left"),
  subset_size = 5,
  ci_level = 0.9,
  precision = 0.001,
  step_fraction = 0.001,
  verbosity = 0,
  round_digits = 3,
  init_tau0 = NULL,
  tau_grid = NULL,
  pvalue_ci = TRUE
)
}
\arguments{
\item{data}{Dataframe of data. This MUST contain the variables refered to in formula.}

\item{coin_test}{Function to calculate test statistic.}

\item{formula}{Formula of outcome on treatment potentiall blocked by some
stratifying factor, e.g., Y ~ Z | D.  Treatment and strat var must be
factors. It is assumed these variables are all inside the passed dataframe.}

\item{tr_level}{Name of tx level to consider "treatment"  Other levels are
control., Default: NULL (and will take first of factor).}

\item{alternative}{What direction (max or minimal effects), Default: 'less'}

\item{distribution}{What kind of ref distribution to calculate (asymp or
approximate), Default: 'asymptotic'}

\item{ci_level}{Level of confidence, Default: 0.9}

\item{step_fraction}{Resolution of grid search, Default: 0.001 (i.e., each
step is this portion of range of Y).}

\item{precision}{How fine should the final grid be for finding the threshold
tau0 for the CI?}

\item{verbosity}{Print out stuff as we go, Default: 1}

\item{round_digits}{Round grid search to this resolution, Default: 3}

\item{init_tau0}{Starting point for search, Default: NULL}

\item{tau_grid}{List of tau values to test in addition to the search.}

\item{pvalue_ci}{Return a confidence interval for the calculated pvalue of
the crossing tau0 (if using approximate MC simulation).}

\item{...}{Passed to the coin_test function}
}
\value{
Confidence interval (one-sided)
}
\description{
Find confidence interval with test inversion

Search a range of values and return a confidence interval for maximal
(minimal) treatment impact.

If there is a stratification variable, permute within strata/blocks defined
by this categorical variable in data, Default: NULL
}
\seealso{
\code{\link[coin]{pvalue}}
}
