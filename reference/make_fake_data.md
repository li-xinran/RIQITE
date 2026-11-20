# Make fake dataset for testing and illustration

These are all fixed set-size and -outcome datasets. The pretend
randomization is perfectly balanced, so we have two units of each type
with one in treatment and one in control.

## Usage

``` r
make_fake_data(scenario = "D")
```

## Arguments

- scenario:

  "A", "B", "C", or "D"; different kinds of data structure.

## Value

Dataset with Y0, Y1, Yobs, and Z (treatment assignment).

## Details

The scenarios are:

A: Long right tail, but also some very low treatment units. Clear
evidence of harm to some units and strong positive effects for others
(the tx units are in the tails, co units are interior)

B: Treatment has lowest values, control has highest. It could be a
simple shift down (constant negative tx effect visually plausible).

C: Treatment helps in general, and helps some a lot.

D: Treatment helps those at the low end, but hurts those at the high
end. So testing for positive tx (alt = greater) on left tail should
reject. Testing for negative tx on right tail(alt = less) should also
reject.
