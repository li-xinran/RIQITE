% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_power.R
\name{calc_power}
\alias{calc_power}
\title{Calculate power for the noninferiority null}
\usage{
calc_power(
  n,
  Y0_distribution = rnorm,
  tx_function = tx_function_factory("constant"),
  p_tx = 0.5,
  R = 100,
  iter_per_set = 10,
  percentile = 1,
  alpha = 0.05,
  c = 0,
  alternative = "greater",
  method.list = list(name = "Stephenson", s = 10),
  score = NULL,
  stat.null = NULL,
  nperm = 1000,
  calc_ICC = FALSE
)
}
\arguments{
\item{n}{Sample size to test}

\item{Y0_distribution}{Either a function to generate a set of Y0
values, or a list of discrete values to sample from (with
replacement).  Function takes single parameter of sample size.
Default is normal.}

\item{tx_function}{Function that takes Y0 and returns a set of
treatment effects.  Default is tx_constant.}

\item{iter_per_set}{Number of permutations per finite dataset.}

\item{calc_ICC}{Calculate how much}
}
\description{
Given a specified DGP for control outcomes and treatment impacts,
calculate power to detect effects.  Calculates superpopulation
power via simulation, where we simulate series of datasets of size
n, randomized into treatment and control, and then test for
detection of effects. Function returns proportion of times effects
are detected.
}
\details{
Each individual dataset is tested iter_per_set times, giving very
rough power calculations for individual sets that are then
aggregated; this provides the ability to do variance decomposition
to see how sensitive power is to individual dataset characteristics
(for datasets of the same family as defined by the DGP) vs.
randomization imbalance.
}
