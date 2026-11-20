# fraction of null draws at least as extreme as the observed statistic

Rank-score statistics are sums of rational scores, and floating point
addition is not associative: two sums that are mathematically tied
(summing the same multiset of scores in a different order) can differ in
their last bit. A strict `>=` comparison then drops a null draw that
should count as tied, undercounting the p-value (see
li-xinran/RIQITE#1). This treats null draws within a small relative
tolerance of `stat.obs` as tied.

## Usage

``` r
pval_from_null(stat.null, stat.obs)
```

## Arguments

- stat.null:

  Vector of the test statistic under the null (e.g. from
  [`null_dist()`](https://li-xinran.github.io/RIQITE/reference/null_dist.md)).

- stat.obs:

  The observed (or hypothesized minimum) value of the test statistic to
  compare the null draws against.
