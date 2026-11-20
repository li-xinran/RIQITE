# Result object for results of ci_quantile call

The ciresult object is an S3 class that holds the results from `ci()`.

It has several methods that pull different information from this object,
and some printing methods for getting nicely formatted results.

Calculate some summary statistics on the result object.

## Usage

``` r
arguments(x)

is.ciresult(x)

# S3 method for class 'ciresult'
head(x, n = 6L, ...)

# S3 method for class 'ciresult'
print(x, ...)

# S3 method for class 'ciresult'
summary(object, c = 0, k = NULL, ...)
```

## Arguments

- x:

  a ciresult object (except for is.ciresult, where it is a generic
  object to check).

- n:

  Number of rows to take

- ...:

  extra options passed to print.ciresult

- object:

  object to summarize.

- c:

  Threshold for testing.

- k:

  Quantiles to report on.

## Value

The parameters of the call

is.ciresult: TRUE if object is a ciresult object.

Data.frame of top rows

print: No return value; prints results.

summary: No return value; prints results.
