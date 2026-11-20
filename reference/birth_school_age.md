# Birth-month and school-entry-age instrumental variable illustration

A reconstructed individual-level dataset used to illustrate testing the
monotonicity of an instrument. This is NOT a set of real individual
records: it is generated to match the aggregate proportions reported in
Black et al. (2011), Table 1 (cross-tabulating month of birth and
school-entry timing) and Table 3 (the "discontinuity subsample" of
104,023 subjects born in January or December, rounded here to 104,000
and split evenly between the two birth months, following Fiorini and
Stevens (2014)).

Subjects born in December tend to enter school on-time or late, while
subjects born in January tend to enter early or on-time, making month of
birth a natural candidate instrument for school-entry age – see
`demo_code/RIllustration20221028.Rmd` for the full worked example
(testing the bounded null and estimating the range of individual effects
of birth month on school-entry age).

## Usage

``` r
birth_school_age
```

## Format

A data frame with 104,000 rows and 3 variables:

- `BirthMonth`:

  character - "Dec" or "Jan".

- `EntryTiming`:

  character - "Early", "On-Time", or "Late".

- `EntryAge`:

  double - school-entry age, in years, implied by `BirthMonth` and
  `EntryTiming`.
