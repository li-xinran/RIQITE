# Explore power across a range of s values

Conduct power simulation against a specific data generating process to
see how power changes as a function of s in that scenario.

## Usage

``` r
explore_stephenson_s(
  s = c(2, 5, 10, 30),
  n,
  Y0_distribution = rnorm,
  tx_function = "constant",
  p_tx = 0.5,
  R = 100,
  iter_per_set = 10,
  percentile = 1,
  alpha = 0.05,
  c = 0,
  k.vec = NULL,
  quantile_n_CI = NA,
  targeted_power = TRUE,
  alternative = "greater",
  nperm = 1000,
  parallel = FALSE,
  n_workers = NULL,
  calc_ICC = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- s:

  Numeric list of s values to test using the Stephenson rank.

- p_tx:

  Proportion of units to be treated.

- R:

  Number of simulation iterations to do.

- percentile:

  For targeted power, what percentile to target.

- alpha:

  Significance level of the test.

- c:

  Threshold bound for treatment impacts for hypothesis testing.

- k.vec:

  Which quantiles to examine (if targeted_power = FALSE).

- quantile_n_CI:

  What quantile of the distribution of lower confidence limits for the
  number of significant units to calculate. NA (the default) means
  calculate mean instead of median. 0.5 would report the median size of
  the confidence interval for number of significant units.

- targeted_power:

  TRUE means focus on passed percentile and calculate power for that
  specific quantile. FALSE means calculate simultenaous confidence
  intervals for all quantiles, leaving power caculation to
  post-simulation exploration.

- alternative:

  A character takes value "greater" and "less", indicating the direction
  of alternative hypotheses.
