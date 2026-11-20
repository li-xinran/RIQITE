# Calculate power or properties of quantile confidence intervals for a finite-sample dataset.

Given a fixed schedule of potential outcomes and a proportion to assign
to treatment, calculate (via simulation) either the power to reject a
null hypothesis or the distribution of all the confidence intervals.

## Usage

``` r
calc_power_finite(
  Y0,
  tau,
  p_tx,
  R = 100,
  percentile = 1,
  alpha = 0.05,
  c = 0,
  alternative = "greater",
  method.list = list(name = "Stephenson", s = 10),
  score = NULL,
  stat.null = NULL,
  use_pval_bound = FALSE,
  nperm = 1000,
  k.vec = NULL,
  quantile_n_CI = NA,
  targeted_power = TRUE,
  summarise = TRUE
)
```

## Arguments

- R:

  Number of permutation replicates for power calculation.

- percentile:

  Which quantile is of interest (default is 1, the max).

- alpha:

  Alpha for testing (default of 0.05)

- c:

  Threshold for null (see pval_quantile)

- alternative:

  A character takes value "greater", "less" and "two.sided", indicating
  the alternative hypothesis.

- method.list:

  A list specifies the choice of the rank sum test statistic. For
  example, list(name="Wilcoxon") means the Wilcoxon rank sum statistic,
  and list(name = "Stephenson", s = 10) means the Stephenson rank sum
  statistic with parameter s=10.

- score:

  An \\n\\ dimensional transformed ranks, i.e., (phi(1), phi(2), ...,
  phi(n)), where phi() denotes the rank transformation function.

- stat.null:

  An vector whose empirical distribution approximates the randomization
  distribution of the rank sum statistic.

- nperm:

  A positive integer representing the number of permutations for
  approximating the randomization distribution of the rank sum
  statistic.

- quantile_n_CI:

  What quantile of the distribution of number of significant units
  should be reported?

- If:

  TRUE, calculates the overall bound on impacts using
  [`pval_bound()`](https://li-xinran.github.io/RIQITE/reference/pval_bound.md)
  (which allows for the Diff in Means test statistic and other stats
  that do not need to be distribution free). Otherwise uses
  [`pval_quantile()`](https://li-xinran.github.io/RIQITE/reference/pval_quantile.md).

## Value

If targeted_power=TRUE, return a small tibble with statistics of the
power calculation (true effect being tested, power to detect, sample
size, etc.). sd_Y0, sd_Y1, sd_tau are the standard deviations of the
potential outcomes and individual treatment impacts, and r is the
correlation of Y0 and the treatment effects. Otherwise, if summarise =
FALSE, return a matrix with each column being a simulation run and each
row being a quantile, with the listed confidence interval limit for each
entry. Otherwise a tibble with a row for each quantile specified in
k.vec. The "n" column is a repeated value of the median number of
quantiles found significant across simulations.#'

## Details

If targeted_power = FALSE, this method will estimate (via simulation)
the distributions all simultaneous confidence intervals that can then be
post processed to determine their variability, and also the CI for the
number of positive units, etc.

Otherwise, it will, by simulation, calculate the chance of rejecting a
null by repeatedly rerandomizing the data and conducting a hypothesis
test on the results, via using
[`pval_quantile()`](https://li-xinran.github.io/RIQITE/reference/pval_quantile.md)
(or
[`pval_bound()`](https://li-xinran.github.io/RIQITE/reference/pval_bound.md)
if `use_pval_bound` is TRUE)
