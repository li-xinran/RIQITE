# Summarize set of power simulation results

Summarize the CI_precision simulation to compare the performance of
various Stephenson rank statistics.

## Usage

``` r
summary_power_results(
  result,
  alternative = "greater",
  c = 0,
  quantile_n_CI = NA
)
```

## Arguments

- result:

  A data frame with one row for each quantile in the result object
  passed (e.g., the result from
  [`calc_power_finite()`](https://li-xinran.github.io/RIQITE/reference/calc_power_finite.md)
  with `targeted_power=FALSE`).

- alternative:

  A character takes value "greater" and "less", indicating the direction
  of alternative hypotheses.

- c:

  A number specifying the cutoff.

- quantile_n_CI:

  What quantile of the distribution of lower confidence limits for the
  number of significant units to calculate. NA (the default) means
  calculate mean instead of median. 0.5 would report the median size of
  the confidence interval for number of significant units.

## Details

Calculates three types of precision. First, we summarize the
distribution of the confidence limit for the number of units with
effects passing the threshold cutoff 'c' (default of c=0). We also, for
each rank, calculate the median confidence limits for the individual
treatment effect for that rank. Finally, we calculate the power for
testing whether the individual treatment effect passes the threshold
cutoff for each rank.
