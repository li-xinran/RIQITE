# Randomization inference for quantiles of individual treatment effects

Obtain one-sided or two-sided confidence intervals for all quantiles of
individual treatment effects

## Usage

``` r
ci_quantile(
  Z,
  Y,
  k.vec = NULL,
  alternative = "greater",
  method.list = list(name = "Stephenson", s = 10),
  score = NULL,
  stat.null = NULL,
  nperm = 10^6,
  Z.perm = NULL,
  alpha = 0.05,
  tol = 10^(-3),
  switch = TRUE
)
```

## Arguments

- Z:

  An \\n\\ dimensional treatment assignment vector.

- Y:

  An \\n\\ dimensional observed outcome vector.

- k.vec:

  A vector that specifies the quantiles of individual effects under
  investigation. If it equals NULL, then we consider all quantiles of
  individual effects.

- alternative:

  A character takes value "greater", "less" and "two.sided", indicating
  whether the confidence intervals are one-sided with lower or upper
  confidence limits or two-sided.

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
  distribution of the rank sum statistic. If NULL, this will be
  calculated internally.

- nperm:

  A positive integer representing the number of permutations for
  approximating the randomization distribution of the rank sum
  statistic.

- Z.perm:

  A \\n \times nperm\\ matrix that specifies the permutated assignments
  for approximating the null distribution of the test statistic.

- alpha:

  A numerical object, where 1-alpha indicates the confidence level.

- tol:

  A numerical object specifying the precision of the obtained confidence
  intervals. For example, if tol = 10^(-3), then the confidence limits
  are precise up to 3 digits.

- switch:

  A logical object indicating whether performing treatment label
  switching to make the treated group have larger size.

## Value

A dataframe with each row corresponding to the lower and upper
confidence limits of the given tau(k). Column names of "lower" and
"upper".
