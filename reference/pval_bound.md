# Randomization test for bounded null hypotheses

Obtain the p-value for testing the bounded null hypothesis H0: tau \<= c
or H0: tau \>= c

## Usage

``` r
pval_bound(
  Z,
  Y,
  c = 0,
  method.list = list(name = "DIM"),
  alternative = "greater",
  Z.perm = NULL,
  nperm = 10^5,
  impute = "control"
)
```

## Arguments

- Z:

  An \\n\\ dimensional treatment assignment vector.

- Y:

  An \\n\\ dimensional observed outcome vector.

- c:

  A scalar or vector specifying the bounded null hypothesis.

- method.list:

  A list specifies the choice of the rank sum test statistic. For
  example, list(name = "DIM") means the difference-in-means statistic,
  list(name="Wilcoxon") means the Wilcoxon rank sum statistic, and
  list(name = "Stephenson", s = 10) means the Stephenson rank sum
  statistic with parameter s=10.

- alternative:

  A character takes value "greater" and "less" indicating the direction
  of the alternative hypothesis.

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

## Value

The p-value for testing the specified null hypothesis of interest.
