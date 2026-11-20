# Randomization test for quantiles of individual treatment effects

Obtain the p-value for testing the null hypothesis H0: \\\tau\_{(k)}
\leq c\\, or \\H0: \tau\_{(k)} \geq c\\, or \\H0: \tau\_{(k)} = c\\,
where \\\tau\_{(k)}\\ denotes individual treatment effect at rank k.

## Usage

``` r
pval_quantile(
  Z,
  Y,
  k,
  c,
  alternative = "greater",
  method.list = list(name = "Stephenson", s = 10),
  score = NULL,
  stat.null = NULL,
  nperm = 10^6,
  Z.perm = NULL,
  switch = TRUE
)
```

## Arguments

- Z:

  An \\n\\ dimensional treatment assignment vector.

- Y:

  An \\n\\ dimensional observed outcome vector.

- k:

  An integer between 1 and n specifying which quantile of individual
  effect is of interest.

- c:

  A numerical object specifying the threshold for the null hypothesis.

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

- Z.perm:

  A \\n \times nperm\\ matrix that specifies the permutated assignments
  for approximating the null distribution of the test statistic.

- switch:

  A logical object indicating whether performing treatment label
  switching to make the treated group have larger size.

## Value

The p-value for testing the specified null hypothesis of interest.
