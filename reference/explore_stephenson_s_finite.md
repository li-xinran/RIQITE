# Explore power and precision across a range of s values

Conduct power simulation against specific finite-sample dataset to see
how power changes as a function of s.

## Usage

``` r
explore_stephenson_s_finite(
  s = c(2, 5, 10, 30),
  Y0,
  tau,
  p_tx,
  R = 100,
  percentile = 1,
  k.vec = NULL,
  alpha = 0.05,
  c = 0,
  quantile_n_CI = NA,
  alternative = "greater",
  stat.nulls = NULL,
  nperm = 1000,
  targeted_power = TRUE
)
```

## Arguments

- s:

  Numeric list of s values to test using the Stephenson rank.

- Y0:

  Either a name of function to generate Y0s, or an empirical
  distribution represented as a list of values.

- p_tx:

  Proportion of units to be treated.

- R:

  Number of simulation iterations to do.

- percentile:

  For targeted power, what percentile to target.

- k.vec:

  Which quantiles to examine (if targeted_power = FALSE).

- alpha:

  Significance level of the test.

- c:

  Threshold bound for treatment impacts for hypothesis testing.

- quantile_n_CI:

  What quantile of the distribution of lower confidence limits for the
  number of significant units to calculate. NA (the default) means
  calculate mean instead of median. 0.5 would report the median size of
  the confidence interval for number of significant units.

- alternative:

  A character takes value "greater" and "less", indicating the direction
  of alternative hypotheses.

- stat.nulls:

  Optional list of the precomputed null distributions of the statistics.

- targeted_power:

  TRUE means focus on passed percentile and calculate power for that
  specific quantile. FALSE means calculate simultenaous confidence
  intervals for all quantiles, leaving power caculation to
  post-simulation exploration.
