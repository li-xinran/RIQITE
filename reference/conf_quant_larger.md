# one-sided confidence interval (c, Inf) for all quantiles

generate one-sided confidence intervals

## Usage

``` r
conf_quant_larger(
  Z,
  Y,
  k.vec = NULL,
  method.list = NULL,
  score = NULL,
  stat.null = NULL,
  nperm = 10^5,
  Z.perm = NULL,
  alpha = 0.05,
  tol = 10^(-3),
  ind.sort.treat = NULL
)
```

## Arguments

- Z:

  treatment assignment (vector)

- Y:

  outcome (vector)
