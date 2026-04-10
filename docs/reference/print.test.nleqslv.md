# Printing the result of `testnslv`

Print method for `test.nleqslv` objects.

## Usage

``` r
# S3 method for class 'test.nleqslv'
print(x, digits = 4, width.cutoff = 45L, ...)
```

## Arguments

- x:

  a `test.nleqslv` object

- digits:

  specifies the minimum number of significant digits to be printed in
  values.

- width.cutoff:

  integer passed to [`deparse`](https://rdrr.io/r/base/deparse.html)
  which sets the cutoff at which line-breaking is tried.

- ...:

  additional arguments to `print`.

## Value

It returns the object `x` invisibly.

## Details

This is the `print` method for objects inheriting from class
`test.nleqslv`. It prints the call to `testnslv` followed by the
description of the experiment (if the `title` argument was specified in
the call to `testnslv`) and the dataframe containing the results of
`testnslv`.

## Examples

``` r
dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}
xstart <- c(1.5,0.5)
fstart <- dslnex(xstart)
z <- testnslv(xstart,dslnex)
print(z)
#> Call:
#> testnslv(x = xstart, fn = dslnex)
#> 
#> Results:
#>     Method Global termcd Fcnt Jcnt Iter Message     Fnorm
#> 1   Newton  cline      1    8    6    6   Fcrit 1.693e-28
#> 2   Newton  qline      1    8    6    6   Fcrit 1.693e-28
#> 3   Newton  gline      1   12    6    6   Fcrit 6.364e-20
#> 4   Newton pwldog      1    8    6    6   Fcrit 2.629e-23
#> 5   Newton dbldog      1    8    6    6   Fcrit 1.710e-26
#> 6   Newton   hook      1    8    6    6   Fcrit 3.753e-24
#> 7   Newton   none      1    7    7    7   Fcrit 8.470e-21
#> 8  Broyden  cline      1   12    1   10   Fcrit 5.315e-19
#> 9  Broyden  qline      1   12    1   10   Fcrit 5.315e-19
#> 10 Broyden  gline      1   19    1   10   Fcrit 3.919e-17
#> 11 Broyden pwldog      1   13    1   11   Fcrit 4.620e-21
#> 12 Broyden dbldog      1   13    1   11   Fcrit 1.041e-21
#> 13 Broyden   hook      1   13    1   10   Fcrit 1.548e-17
#> 14 Broyden   none      1   13    1   13   Fcrit 8.286e-18
```
