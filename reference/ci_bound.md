# Confidence interval for the maximum or minimum individual effect

Obtain confidence interval for the maximum or minimum individual effect

## Usage

``` r
ci_bound(
  Z,
  Y,
  alternative = "greater",
  method.list = list(name = "DIM"),
  Z.perm = NULL,
  nperm = 10^3,
  impute = "control",
  alpha = 0.05,
  tol = 10^(-3)
)
```

## Arguments

- Z:

  An \\n\\ dimensional treatment assignment vector.

- Y:

  An \\n\\ dimensional observed outcome vector.

- alternative:

  A character takes value "greater" and "less" indicating the direction
  of the alternative hypothesis.

- method.list:

  A list specifies the choice of the rank sum test statistic. For
  example, list(name = "DIM") means the difference-in-means statistic,
  list(name="Wilcoxon") means the Wilcoxon rank sum statistic, and
  list(name = "Stephenson", s = 10) means the Stephenson rank sum
  statistic with parameter s=10.

- Z.perm:

  A \\n \times nperm\\ matrix that specifies the permutated assignments
  for approximating the null distribution of the test statistic.

- nperm:

  A positive integer representing the number of permutations for
  approximating the randomization distribution of the test statistic

- impute:

  A character specifies the form of the test statistic. "control" means
  we use the imputed control potential outcomes in the test statistic,
  and "both" means we use the imputed observed outcome in the test
  statistic.

- alpha:

  A numeric object specifies the level of the confidence interval. In
  particular, the confidence level is \\1-\alpha\\.

- tol:

  A numerical object specifying the precision of the obtained confidence
  intervals. For example, if tol = 10^(-3), then the confidence limits
  are precise up to 3 digits.

## Value

The p-value for testing the specified null hypothesis of interest.
