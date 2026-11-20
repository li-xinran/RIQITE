# Generate randomization dist of the rank score stat

Generate the null distribution of the given test statistic for an
experiment with m treated out of n units.

## Usage

``` r
null_dist(n, m, method.list = NULL, score = NULL, nperm = 10^5, Z.perm = NULL)
```

## Arguments

- score:

  The list of possible scores of the units (under a null, invariant
  under treatment assignment). If null will generate from the passed
  method.list and given n, m.

- Z.perm:

  Matrix of possible treatment assignments. NULL will generate complete
  rand with nperm combinations.
