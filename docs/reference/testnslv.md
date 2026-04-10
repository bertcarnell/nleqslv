# Test different methods for solving with `nleqslv`

The function tests different methods and global strategies for solving a
system of nonlinear equations with `nleqslv`

## Usage

``` r
testnslv(
  x,
  fn,
  jac = NULL,
  ...,
  method = c("Newton", "Broyden"),
  global = c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none"),
  Nrep = 0L,
  title = NULL
)
```

## Arguments

- x:

  A numeric vector with an initial guess of the root.

- fn:

  A function of `x` returning the function values.

- jac:

  A function to return the Jacobian for the `fn` function. For a vector
  valued function `fn` the Jacobian must be a numeric matrix of the
  correct dimensions. For a scalar valued function `fn` the `jac`
  function may return a scalar. If not supplied numerical derivatives
  will be used.

- ...:

  Further arguments to be passed to `fn` and `jac` and
  [`nleqslv`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.md).

- method:

  The methods to use for finding a solution.

- global:

  The global strategies to test. The argument may consist of several
  possibly abbreviated items.

- Nrep:

  Number of repetitions to apply. Default is no repetitions.

- title:

  a description of this experiment.

## Value

`testnslv` returns an object of class `"test.nleqslv"` which is a list
containing the following elements

- `call`:

  the matched call

- `out`:

  a dataframe containing the results with the following columns

  `Method`

  :   method used.

  `Global`

  :   global strategy used.

  `termcd`

  :   termination code of `nleqslv`.

  `Fcnt`

  :   number of function evaluations used by the method and global
      strategy. This excludes function evaluations made when computing a
      numerical Jacobian.

  `Jcnt`

  :   number of Jacobian evaluations.

  `Iter`

  :   number of outer iterations used by the algorithm.

  `Message`

  :   a string describing the termination code in an abbreviated form.

  `Fnorm`

  :   square of the euclidean norm of the vector of function results
      divided by 2.

  `cpusecs`

  :   CPU seconds used by the requested number of repetitions (only
      present when argument `Nrep` is not 0).

- `title`:

  the description if specified

The abbreviated strings are

- `Fcrit`:

  Convergence of function values has been achieved.

- `Xcrit`:

  This means that the relative distance between two consecutive x-values
  is smaller than `xtol`.

- `Stalled`:

  The algorithm cannot find an acceptable new point.

- `Maxiter`:

  Iteration limit `maxit` exceeded.

- `Illcond`:

  Jacobian is too ill-conditioned.

- `Singular`:

  Jacobian is singular.

- `BadJac`:

  Jacobian is unusable.

- `ERROR`:

  `nleqslv` stopped because of a fatal error.

## Details

The function solves the function `fn` with
[`nleqslv`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.md)
for the specified methods and global strategies. When argument `Nrep`
has been set to a number greater than or equal to 1, repetitions of the
solving process are performed and the used CPU time in seconds is
recorded.

If checking a user supplied jacobian is enabled, then `testnslv` will
stop immediately when a possibly incorrect jacobian is detected.

## Warning

Any `nleqslv` error message will be displayed immediately and an error
for the particular combination of method and global strategy will be
recorded in the final dataframe.

## Examples

``` r
dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}
xstart <- c(0.5,0.5)
fstart <- dslnex(xstart)
testnslv(xstart,dslnex)
#> Call:
#> testnslv(x = xstart, fn = dslnex)
#> 
#> Results:
#>     Method Global termcd Fcnt Jcnt Iter Message     Fnorm
#> 1   Newton  cline      1    7    6    6   Fcrit 8.550e-23
#> 2   Newton  qline      1    7    6    6   Fcrit 8.550e-23
#> 3   Newton  gline      1    9    5    5   Fcrit 1.501e-20
#> 4   Newton pwldog      1    7    6    6   Fcrit 4.536e-27
#> 5   Newton dbldog      1    7    6    6   Fcrit 4.536e-27
#> 6   Newton   hook      1    7    6    6   Fcrit 2.450e-23
#> 7   Newton   none      1    8    8    8   Fcrit 6.085e-25
#> 8  Broyden  cline      1   12    1    9   Fcrit 4.297e-17
#> 9  Broyden  qline      1   12    1    9   Fcrit 4.297e-17
#> 10 Broyden  gline      1   14    1   10   Fcrit 4.960e-21
#> 11 Broyden pwldog      1   12    1   10   Fcrit 5.156e-21
#> 12 Broyden dbldog      1   12    1   10   Fcrit 3.958e-20
#> 13 Broyden   hook      1   12    1   10   Fcrit 1.127e-19
#> 14 Broyden   none      1   13    1   13   Fcrit 1.783e-22
# this will encounter an error
xstart <- c(2.0,0.5)
fstart <- dslnex(xstart)
testnslv(xstart,dslnex)
#> Error (method=Newton global=none): non-finite value(s) detected in jacobian (row=2,col=1)
#> Call:
#> testnslv(x = xstart, fn = dslnex)
#> 
#> Results:
#>     Method Global termcd Fcnt Jcnt Iter Message     Fnorm
#> 1   Newton  cline      1   11    7    7   Fcrit 4.992e-19
#> 2   Newton  qline      1   10    7    7   Fcrit 2.809e-20
#> 3   Newton  gline      1   17    7    7   Fcrit 6.311e-29
#> 4   Newton pwldog      1    6    5    5   Fcrit 1.543e-18
#> 5   Newton dbldog      1    6    5    5   Fcrit 1.790e-18
#> 6   Newton   hook      1   11    7    7   Fcrit 3.819e-26
#> 7   Newton   none     NA   NA   NA   NA   ERROR        NA
#> 8  Broyden  cline      1   17    1   11   Fcrit 2.900e-18
#> 9  Broyden  qline      1   18    1   13   Fcrit 1.404e-17
#> 10 Broyden  gline      1   25    1   11   Fcrit 4.798e-19
#> 11 Broyden pwldog      1   12    1   10   Fcrit 6.230e-18
#> 12 Broyden dbldog      1   12    1   10   Fcrit 3.239e-18
#> 13 Broyden   hook      1   16    1   12   Fcrit 5.719e-23
#> 14 Broyden   none      4   20    1   20 Maxiter 1.846e-01
```
