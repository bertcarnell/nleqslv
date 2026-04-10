# Solve a nonlinear equation system with multiple roots from multiple initial estimates

This function solves a system of nonlinear equations with `nleqlsv` for
multiple initial estimates of the roots.

## Usage

``` r
searchZeros(x, fn, digits = 4L, ...)
```

## Arguments

- x:

  A matrix with each row containing an initial guess of the roots.

- fn:

  A function of `x` returning a vector of function values with the same
  length as the vector `x`.

- digits:

  integer passed to [`round`](https://rdrr.io/r/base/Round.html) for
  locating and removing duplicate rounded solutions.

- ...:

  Further arguments to be passed to
  [`nleqslv`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.md),
  `fn` and `jac`.

## Value

If no solutions are found `NULL` is returned. Otherwise a list
containing the following components is returned

- `x`:

  a matrix with each row containing a unique solution (unrounded)

- `xfnorm`:

  a vector of the function criterion associated with each row of the
  solution matrix `x`.

- `fnorm`:

  a vector containing the function criterion for every converged result

- `idxcvg`:

  a vector containing the row indices of the matrix of initial estimates
  for which function value convergence was achieved

- `idxxtol`:

  a vector containing the row indices of the matrix of initial estimates
  for which x-value convergence was achieved

- `idxnocvg`:

  a vector containing the row indices of the matrix of initial estimates
  which lead to an `nleqslv` termination code \> 2

- `idxfatal`:

  a vector containing the row indices of the matrix of initial estimates
  for which a fatal error occurred that made `nleqslv` stop

- `xstart`:

  a matrix of the initial estimates corresponding to the solution matrix

- `cvgstart`:

  a matrix of all initial estimates for which convergence was achieved

## Details

Each row of `x` is a vector of initial estimates for the argument `x` of
`nleqslv`. The function runs `nleqslv` for each row of the matrix `x`.
The first initial value is treated separately and slightly differently
from the other initial estimates. It is used to check if all arguments
in `...` are valid arguments for `nleqslv` and the function to be
solved. This is done by running `nleqslv` with no condition handling. If
an error is then detected an error message is issued and the function
stops. For the remaining initial estimates `nleqslv` is executed
silently. Only solutions for which the `nleqslv` termination code
`tcode` equals `1` are regarded as valid solutions. The rounded
solutions (after removal of duplicates) are used to order the solutions
in increasing order. These rounded solutions are not included in the
return value of the function.

## Examples

``` r
# Dennis Schnabel example 6.5.1 page 149 (two solutions)
set.seed(123)
dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}
xstart <- matrix(runif(50, min=-2, max=2),ncol=2)
ans <- searchZeros(xstart,dslnex, method="Broyden",global="dbldog")
ans
#> $x
#>            [,1]     [,2]
#> [1,] -0.7137474 1.220887
#> [2,]  1.0000000 1.000000
#> 
#> $xfnorm
#> [1] 1.180429e-20 1.799915e-22
#> 
#> $fnorm
#>  [1] 1.180429e-20 8.445681e-20 1.844667e-21 1.799915e-22 4.360810e-20
#>  [6] 4.822898e-19 5.148702e-18 4.267196e-21 8.063623e-18 2.137277e-21
#> [11] 5.381588e-19 5.259062e-23 1.742763e-18
#> 
#> $idxcvg
#>  [1]  1  3  6  7  8  9 10 12 15 17 18 19 25
#> 
#> $idxxtol
#> integer(0)
#> 
#> $idxnocvg
#>  [1]  2  4  5 11 13 14 16 20 21 22 23 24
#> 
#> $idxfatal
#> integer(0)
#> 
#> $xstart
#>            [,1]      [,2]
#> [1,] -0.8496899 0.8341219
#> [2,]  0.1124220 1.6091962
#> 
#> $cvgstart
#>             [,1]       [,2]
#>  [1,] -0.8496899  0.8341219
#>  [2,] -0.3640923  0.3765681
#>  [3,] -1.8177740  1.8520969
#>  [4,]  0.1124220  1.6091962
#>  [5,]  1.5696762  0.7628211
#>  [6,]  0.2057401  1.1818697
#>  [7,] -0.1735411 -1.9015453
#>  [8,] -0.1866634  1.0338382
#>  [9,] -1.5883013 -1.0734969
#> [10,] -1.0156491 -0.3418147
#> [11,] -1.8317619 -0.3451027
#> [12,] -0.6883171 -0.5246182
#> [13,]  0.6228232  1.4313109
#> 

# more complicated example
# R. Baker Kearfott, Some tests of Generalized Bisection,
# ACM Transactions on Methematical Software, Vol. 13, No. 3, 1987, pp 197-220

# A high-degree polynomial system (section 4.3 Problem 12)
# There are 12 real roots (and 126 complex roots to this system!)

hdp <- function(x) {
    f <- numeric(length(x))
    f[1] <- 5 * x[1]^9 - 6 * x[1]^5 * x[2]^2 + x[1] * x[2]^4 + 2 * x[1] * x[3]
    f[2] <- -2 * x[1]^6 * x[2] + 2 * x[1]^2 * x[2]^3 + 2 * x[2] * x[3]
    f[3] <- x[1]^2 + x[2]^2 - 0.265625
    f
}


N <- 40 # at least to find all 12 roots
set.seed(123)
xstart <- matrix(runif(3*N,min=-1,max=1), N, 3)  # N initial guesses, each of length 3
ans <- searchZeros(xstart,hdp, method="Broyden",global="dbldog")
ans$x
#>                [,1]          [,2]          [,3]
#>  [1,] -5.153882e-01  1.399731e-09 -1.244560e-02
#>  [2,] -4.669800e-01 -2.180703e-01 -2.288440e-10
#>  [3,] -4.669801e-01  2.180702e-01  7.797307e-10
#>  [4,] -2.798547e-01 -4.327890e-01 -1.418919e-02
#>  [5,] -2.798547e-01  4.327890e-01 -1.418919e-02
#>  [6,] -1.689587e-08 -5.153882e-01 -1.484649e-10
#>  [7,]  3.838786e-09  5.153882e-01 -3.274405e-11
#>  [8,]  2.798547e-01 -4.327890e-01 -1.418919e-02
#>  [9,]  2.798547e-01  4.327890e-01 -1.418919e-02
#> [10,]  4.669800e-01 -2.180703e-01 -2.864442e-10
#> [11,]  4.669800e-01  2.180704e-01  6.597637e-10
#> [12,]  5.153882e-01 -3.239027e-09 -1.244560e-02
```
