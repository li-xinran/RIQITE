# Calculate power for the noninferiority null

Given a specified DGP for control outcomes and treatment impacts,
calculate power to detect effects. Calculates superpopulation power via
simulation, where we simulate series of datasets of size n, randomized
into treatment and control, and then test for detection of effects.
Function returns proportion of times effects are detected.

## Usage

``` r
calc_power(
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
  alternative = "greater",
  method.list = list(name = "Stephenson", s = 10),
  targeted_power = TRUE,
  use_pval_bound = FALSE,
  score = NULL,
  stat.null = NULL,
  nperm = 1000,
  calc_ICC = FALSE,
  quantile_n_CI = NA,
  summarise = TRUE,
  ...
)
```

## Arguments

- n:

  Sample size to test

- Y0_distribution:

  Either a function to generate a set of Y0 values, or a list of
  discrete values to bootstrap sample from (with replacement). Function
  takes single parameter of sample size. Default is normal.

- tx_function:

  Either a string name of an implemented distribution (rexp, rnorm,
  linear, etc.) or a function that takes Y0 and returns a set of
  treatment effects. Default is a constant tx effect. tx_constant.

- R:

  Number of permutation replicates for power calculation.

- iter_per_set:

  Number of permutations per finite dataset.

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

- calc_ICC:

  Calculate how much power varies across finite datasets.

- quantile_n_CI:

  What quantile of the distribution of number of significant units
  should be reported?

- ...:

  Extra parameters to generate treatment impacts with.

## Details

Each individual dataset is tested iter_per_set times, giving very rough
power calculations for individual sets that are then aggregated; this
provides the ability to do variance decomposition to see how sensitive
power is to individual dataset characteristics (for datasets of the same
family as defined by the DGP) vs. randomization imbalance.
