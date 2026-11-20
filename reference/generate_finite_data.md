# Generate data for power simulation

Given a desired sample size and a data generating mechanism, generate a
dataset of Y0s and individual treatment effects.

## Usage

``` r
generate_finite_data(n, Y0_distribution = rnorm, tx_function = "constant", ...)
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

- ...:

  Extra parameters to generate treatment impacts with.

## Value

Tibble (dataframe) with two columns, one of Y0s and one (tau) of
treatment effects.
