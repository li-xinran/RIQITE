# Teachers of Electric Circuits Professional Development

This dataset describes the baseline and final test scores on a test
measuring the knowledge of electric circuits given to teachers who were
part of a randomized controlled trial of a professional development
program.

The trial had three treatment arms (different approaches to the
professional development) and a control arm.

The gain (post score - pre score) for all teachers is the change in the
percentage point correct on the test of electric circuit knowledge (we
call this "content knowledge").

See Heller, J. I., Daehler, K. R., Wong, N., Shinohara, M., & Miratrix,
L. W. (2012). Differential effects of three professional development
models on teacher knowledge and student achievement in elementary
science. Journal of Research in Science Teaching, 49(3), 333–362. doi:
10.1002/tea.21004

See report on ERIC at https://files.eric.ed.gov/fulltext/ED514193.pdf

## Usage

``` r
electric_teachers
```

## Format

A data frame with 233 rows and 7 variables:

- `Site`:

  character - Site ID (city where teachers were recruited).
  Randomization occured within site.

- `T.ID`:

  double - Teacher ID

- `Tx`:

  character - Treatment group. "D" is control.

- `TxAny`:

  double - Flag for any treatment (main analysis did not find major
  differences of PD approach.)

- `Per.1`:

  double - Baseline knowledge (percent correct)

- `Per.2`:

  double - Endline knowledge (percent correct)

- `gain`:

  double - Total gain from pre to post.

## Details

Not all treatments were given in a site; PD program was delivered to
classes of teachers. Teachers individually randomized to class, but
there could be class (or PD teacher) effects.
