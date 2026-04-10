# Solving systems of nonlinear equations with Broyden or Newton

The function solves a system of nonlinear equations with either a
Broyden or a full Newton method. It provides line search and trust
region global strategies for difficult systems.

## Usage

``` r
nleqslv(
  x,
  fn,
  jac = NULL,
  ...,
  method = c("Broyden", "Newton"),
  global = c("dbldog", "pwldog", "cline", "qline", "gline", "hook", "none"),
  xscalm = c("fixed", "auto"),
  jacobian = FALSE,
  control = list()
)
```

## Arguments

- x:

  A numeric vector with an initial guess of the root of the function.

- fn:

  A function of `x` returning a vector of function values with the same
  length as the vector `x`.

- jac:

  A function to return the Jacobian for the `fn` function. For a vector
  valued function `fn` the Jacobian must be a numeric matrix of the
  correct dimensions. For a scalar valued function `fn` the `jac`
  function may return a scalar. If not supplied numerical derivatives
  will be used.

- ...:

  Further arguments to be passed to `fn` and `jac`.

- method:

  The method to use for finding a solution. See ‘Details’.

- global:

  The global strategy to apply. See ‘Details’.

- xscalm:

  The type of x scaling to use. See ‘Details’.

- jacobian:

  A logical indicating if the estimated (approximate) jacobian in the
  solution should be returned. See ‘Details’.

- control:

  A named list of control parameters. See ‘Details’.

## Value

A list containing components

- x:

  final values for x

- fvec:

  function values

- termcd:

  termination code as integer. The values returned are

  `1`

  :   Function criterion is near zero. Convergence of function values
      has been achieved.

  `2`

  :   x-values within tolerance. This means that the relative distance
      between two consecutive x-values is smaller than `xtol` but that
      the function value criterion is still larger than `ftol`.
      *Function values may not be near zero; therefore the user must
      check if function values are acceptably small.*

  `3`

  :   No better point found. This means that the algorithm has stalled
      and cannot find an acceptable new point. This may or may not
      indicate acceptably small function values.

  `4`

  :   Iteration limit `maxit` exceeded.

  `5`

  :   Jacobian is too ill-conditioned.

  `6`

  :   Jacobian is singular.

  `7`

  :   Jacobian is unusable.

  `-10`

  :   User supplied Jacobian is most likely incorrect.

- message:

  a string describing the termination code

- scalex:

  a vector containing the scaling factors, which will be the final
  values when automatic scaling was selected

- njcnt:

  number of Jacobian evaluations

- nfcnt:

  number of function evaluations, excluding those required for
  calculating a Jacobian and excluding the initial function evaluation
  (at iteration 0)

- iter:

  number of outer iterations used by the algorithm. This excludes the
  initial iteration. The number of backtracks can be calculated as the
  difference between the `nfcnt` and `iter` components.

- jac:

  the final Jacobian or the Broyden approximation if `jacobian` was set
  to `TRUE`. If no iterations were executed, as may happen when the
  initial guess are sufficiently close the solution, there is no Broyden
  approximation and the returned matrix will always be the actual
  Jacobian. If the matrix is singular or too ill-conditioned the
  returned matrix is of no value.

## Details

The algorithms implemented in `nleqslv` are based on Dennis and Schnabel
(1996).

### Methods

Method `Broyden` starts with a computed Jacobian of the function and
then updates this Jacobian after each successful iteration using the
so-called Broyden update. This method often shows super linear
convergence towards a solution. When `nleqslv` determines that it cannot
continue with the current Broyden matrix it will compute a new Jacobian.

Method `Newton` calculates a Jacobian of the function `fn` at each
iteration. Close to a solution this method will usually show quadratic
convergence.

Both methods apply a so-called (backtracking) global strategy to find a
better (more acceptable) iterate. The function criterion used by the
algorithm is half of the sum of squares of the function values and
“acceptable” means sufficient decrease of the current function criterion
value compared to that of the previous iteration. A comprehensive
discussion of these issues can be found in Dennis and Schnabel (1996).
Both methods apply an unpivoted QR-decomposition to the Jacobian as
implemented in Lapack. The Broyden method applies a rank-1 update to the
Jacobian at the end of each iteration and is based on a simplified and
modernized version of the algorithm described in Reichel and Gragg
(1990).

### Global strategies

When applying a full Newton or Broyden step does not yield a
sufficiently smaller function criterion value `nleqslv` will attempt to
decrease the steplength using one of several so-called global
strategies.

The `global` argument indicates which global strategy to use or to use
no global strategy

- `cline`:

  a cubic line search

- `qline`:

  a quadratic line search

- `gline`:

  a geometric line search

- `dbldog`:

  a trust region method using the double dogleg method as described in
  Dennis and Schnabel (1996)

- `pwldog`:

  a trust region method using the Powell dogleg method as developed by
  Powell (1970).

- `hook`:

  a trust region method described by Dennis and Schnabel (1996) as *The
  locally constrained optimal (“hook”) step*. It is equivalent to a
  Levenberg-Marquardt algorithm as described in MoréMore (1978) and
  Nocedal and Wright (2006).

- `none`:

  Only a pure local Newton or Broyden iteration is used. The maximum
  stepsize (see below) is taken into account. The default maximum number
  of iterations (see below) is set to 20.

The double dogleg method is the default global strategy employed by this
package.

Which global strategy to use in a particular situation is a matter of
trial and error. When one of the trust region methods fails, one of the
line search strategies should be tried. Sometimes a trust region will
work and sometimes a line search method; neither has a clear advantage
but in many cases the double dogleg method works quite well.

When the function to be solved returns non-finite function values for a
parameter vector `x` and the algorithm is *not* evaluating a numerical
Jacobian, then any non-finite values will be replaced by a large number
forcing the algorithm to backtrack, i.e. decrease the line search factor
or decrease the trust region radius.

### Scaling

The elements of vector `x` may be scaled during the search for a zero of
`fn`. The `xscalm` argument provides two possibilities for scaling

- `fixed`:

  the scaling factors are set to the values supplied in the `control`
  argument and remain unchanged during the iterations. The scaling
  factor of any element of `x` should be set to the inverse of the
  typical value of that element of `x`, ensuring that all elements of
  `x` are approximately equal in size.

- `auto`:

  the scaling factors are calculated from the euclidean norms of the
  columns of the Jacobian matrix. When a new Jacobian is computed, the
  scaling values will be set to the euclidean norm of the corresponding
  column if that is larger than the current scaling value. Thus the
  scaling values will not decrease during the iteration. This is the
  method described in MoréMore (1978). Usually manual scaling is
  preferable.

### Jacobian

When evaluating a numerical Jacobian, an error message will be issued on
detecting non-finite function values. An error message will also be
issued when a user supplied jacobian contains non-finite entries.

When the `jacobian` argument is set to `TRUE` the final Jacobian or
Broyden matrix will be returned in the return list. The default value is
`FALSE`; i.e. to not return the final matrix. There is no guarantee that
the final Broyden matrix resembles the actual Jacobian.

The package can cope with a singular or ill-conditioned Jacobian if
needed by setting the `allowSingular` component of the `control`
argument. The method used is described in Dennis and Schnabel (1996); it
is equivalent to a Levenberg-Marquardt type adjustment with a small
damping factor. *There is no guarantee that this method will be
successful.* Warning: *`nleqslv` may report spurious convergence in this
case.*

By default `nleqslv` returns an error if a Jacobian becomes singular or
very ill-conditioned. A Jacobian is considered to be very
ill-conditioned when the estimated inverse condition is less than or
equal to a specified tolerance with a default value equal to
\\10^{-12}\\; this can be changed and made smaller with the `cndtol`
item of the `control` argument. *There is no guarantee that any change
will be effective.*

### Control options

The `control` argument is a named list that can supply any of the
following components:

- `xtol`:

  The relative steplength tolerance. When the relative steplength of all
  scaled x values is smaller than this value convergence is declared.
  The default value is \\10^{-8}\\.

- `ftol`:

  The function value tolerance. Convergence is declared when the largest
  absolute function value is smaller than `ftol`. The default value is
  \\10^{-8}\\.

- `btol`:

  The backtracking tolerance. When the relative steplength in a
  backtracking step to find an acceptable point is smaller than the
  backtracking tolerance, the backtracking is terminated. In the
  `Broyden` method a new Jacobian will be calculated if the Jacobian is
  outdated. The default value is \\10^{-3}\\.

- `cndtol`:

  The tolerance of the test for ill conditioning of the Jacobian or
  Broyden approximation. If less than the machine precision it will be
  silently set to the machine precision. When the estimated inverse
  condition of the (approximated) Jacobian matrix is less than or equal
  to the value of `cndtol` the matrix is deemed to be ill-conditioned,
  in which case an error will be reported if the `allowSingular`
  component is set to `FALSE`. The default value is \\10^{-12}\\.

- `sigma`:

  Reduction factor for the geometric line search. The default value is
  0.5.

- `scalex`:

  a vector of scaling values for the parameters. The inverse of a scale
  value is an indication of the size of a parameter. The default value
  is 1.0 for all scale values.

- `maxit`:

  The maximum number of major iterations. The default value is 150 if a
  global strategy has been specified. If no global strategy has been
  specified the default is 20.

- `trace`:

  Non-negative integer. A value of 1 will give a detailed report of the
  progress of the iteration. For a description see
  [`Iteration-report`](https://bertcarnell.github.io/nleqslv/reference/nleqslv-iterationreport.md).

- `chkjac`:

  A logical value indicating whether to check a user supplied Jacobian,
  if supplied. The default value is `FALSE`. The first 10 errors are
  printed. The code for this check is derived from the code in Bouaricha
  and Schnabel (1997).

- `delta`:

  Initial (scaled) trust region radius. A value of \\-1.0\\ or
  `"cauchy"` is replaced by the length of the Cauchy step in the initial
  point. A value of \\-2.0\\ or `"newton"` is replaced by the length of
  the Newton step in the initial point. Any numeric value less than or
  equal to 0 and not equal to \\-2.0\\, will be replaced by \\-1.0\\;
  the algorithm will then start with the length of the Cauchy step in
  the initial point. If it is numeric and positive it will be set to the
  smaller of the value supplied or the maximum stepsize. If it is not
  numeric and not one of the permitted character strings then an error
  message will be issued. The default is \\-2.0\\.

- `stepmax`:

  Maximum scaled stepsize. If this is negative then the maximum stepsize
  is set to the largest positive representable number. The default is
  \\-1.0\\, so there is no default maximum stepsize.

- `dsub`:

  Number of non zero subdiagonals of a banded Jacobian. The default is
  to assume that the Jacobian is *not* banded. Must be specified if
  `dsuper` has been specified and must be larger than zero when `dsuper`
  is zero.

- `dsuper`:

  Number of non zero super diagonals of a banded Jacobian. The default
  is to assume that the Jacobian is *not* banded. Must be specified if
  `dsub` has been specified and must be larger than zero when `dsub` is
  zero.

- `allowSingular`:

  A logical value indicating if a small correction to the Jacobian when
  it is singular or too ill-conditioned is allowed. If the correction is
  less than `100*.Machine$double.eps` the correction cannot be applied
  and an unusable Jacobian will be reported. The method used is similar
  to a Levenberg-Marquardt correction and is explained in Dennis and
  Schnabel (1996) on page 151. It may be necessary to choose a higher
  value for `cndtol` to enforce the correction. The default value is
  `FALSE`.

## Warning

You cannot use this function recursively. Thus function `fn` should not
in its turn call `nleqslv`.

## References

Bouaricha, A. and Schnabel, R.B. (1997), Algorithm 768: TENSOLVE: A
Software Package for Solving Systems of Nonlinear Equations and
Nonlinear Least-squares Problems Using Tensor Methods, *Transactions on
Mathematical Software*, **23**, 2, pp. 174–195.

Dennis, J.E. Jr and Schnabel, R.B. (1996), *Numerical Methods for
Unconstrained Optimization and Nonlinear Equations*, Siam.

MoréMore, J.J. (1978), The Levenberg-Marquardt Algorithm, Implementation
and Theory, In *Numerical Analysis*, G.A. Watson (Ed.), Lecture Notes in
Mathematics 630, Springer-Verlag, pp. 105–116.

Golub, G.H and C.F. Van Loan (1996), Matrix Computations (3rd edition),
The John Hopkins University Press.

Higham, N.J. (2002), Accuracy and Stability of Numerical Algorithms, 2nd
ed., SIAM, pp. 10–11.

Nocedal, J. and Wright, S.J. (2006), *Numerical Optimization*, Springer.

Powell, M.J.D. (1970), A hybrid method for nonlinear algebraic
equations, In *Numerical Methods for Nonlinear Algebraic Equations*, P.
Rabinowitz (Ed.), Gordon & Breach.

Powell, M.J.D. (1970), A Fortran subroutine for solving systems
nonlinear equations, In *Numerical Methods for Nonlinear Algebraic
Equations*, P. Rabinowitz (Ed.), Gordon & Breach.

decompositions,

Reichel, L. and W.B. Gragg (1990), Algorithm 686: FORTRAN subroutines
for updating the QR decomposition, *ACM Trans. Math. Softw.*, **16**, 4,
pp. 369–377.

## See also

If this function cannot solve the supplied function then it is a good
idea to try the function
[testnslv](https://bertcarnell.github.io/nleqslv/reference/testnslv.md)
in this package. For detecting multiple solutions see
[searchZeros](https://bertcarnell.github.io/nleqslv/reference/searchZeros.md).

## Examples

``` r
# Dennis Schnabel example 6.5.1 page 149
dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}

jacdsln <- function(x) {
    n <- length(x)
    Df <- matrix(numeric(n*n),n,n)
    Df[1,1] <- 2*x[1]
    Df[1,2] <- 2*x[2]
    Df[2,1] <- exp(x[1]-1)
    Df[2,2] <- 3*x[2]^2

    Df
}

BADjacdsln <- function(x) {
    n <- length(x)
    Df <- matrix(numeric(n*n),n,n)
    Df[1,1] <- 4*x[1]
    Df[1,2] <- 2*x[2]
    Df[2,1] <- exp(x[1]-1)
    Df[2,2] <- 5*x[2]^2

    Df
}

xstart <- c(2,0.5)
fstart <- dslnex(xstart)
xstart
#> [1] 2.0 0.5
fstart
#> [1] 2.2500000 0.8432818

# a solution is c(1,1)

nleqslv(xstart, dslnex, control=list(btol=.01))
#> $x
#> [1] 1 1
#> 
#> $fvec
#> [1] 1.499777e-09 2.056389e-09
#> 
#> $termcd
#> [1] 1
#> 
#> $message
#> [1] "Function criterion near zero"
#> 
#> $scalex
#> [1] 1 1
#> 
#> $nfcnt
#> [1] 12
#> 
#> $njcnt
#> [1] 1
#> 
#> $iter
#> [1] 10
#> 

# Cauchy start
nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="cauchy"))
#>   Algorithm parameters
#>   --------------------
#>   Method: Broyden  Global strategy: double dogleg (initial trust region = -1)
#>   Maximum stepsize = 1.79769e+308
#>   Scaling: fixed
#>   ftol = 1e-08 xtol = 1e-08 btol = 0.01 cndtol = 1e-12
#> 
#>   Iteration report
#>   ----------------
#>   Iter         Jac     Lambda      Eta     Dlt0     Dltn         Fnorm   Largest |f|
#>      0                                                    2.886812e+00  2.250000e+00
#>      1  N(9.6e-03) C            0.9544   0.4671   0.9343* 1.699715e-01  5.421673e-01
#>      1             W   0.0833   0.9544   0.9343   0.4671  1.699715e-01  5.421673e-01
#>      2  B(1.1e-02) W   0.1154   0.4851   0.4671   0.4671  1.277667e-01  5.043571e-01
#>      3  B(7.3e-02) W   0.7879   0.7289   0.4671   0.0759  5.067893e-01  7.973542e-01
#>      3             C            0.7289   0.0759   0.1519  5.440250e-02  2.726084e-01
#>      4  B(8.3e-02) W   0.5307   0.3271   0.1519   0.3037  3.576547e-02  2.657553e-01
#>      5  B(1.8e-01) N            0.6674   0.2191   0.4383  6.566182e-03  8.555110e-02
#>      6  B(1.8e-01) N            0.9801   0.0376   0.0752  4.921645e-04  3.094104e-02
#>      7  B(1.9e-01) N            0.7981   0.0157   0.0313  4.960629e-06  2.826064e-03
#>      8  B(1.6e-01) N            0.3942   0.0029   0.0058  1.545503e-08  1.757498e-04
#>      9  B(1.5e-01) N            0.6536   0.0001   0.0003  2.968676e-11  5.983765e-06
#>     10  B(1.5e-01) N            0.4730   0.0000   0.0000  4.741792e-14  2.198380e-07
#>     11  B(1.5e-01) N            0.9787   0.0000   0.0000  6.451792e-19  8.118586e-10
#> $x
#> [1] 1 1
#> 
#> $fvec
#> [1]  8.118586e-10 -7.945087e-10
#> 
#> $termcd
#> [1] 1
#> 
#> $message
#> [1] "Function criterion near zero"
#> 
#> $scalex
#> [1] 1 1
#> 
#> $nfcnt
#> [1] 13
#> 
#> $njcnt
#> [1] 1
#> 
#> $iter
#> [1] 11
#> 

# Newton start
nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))
#>   Algorithm parameters
#>   --------------------
#>   Method: Broyden  Global strategy: double dogleg (initial trust region = -2)
#>   Maximum stepsize = 1.79769e+308
#>   Scaling: fixed
#>   ftol = 1e-08 xtol = 1e-08 btol = 0.01 cndtol = 1e-12
#> 
#>   Iteration report
#>   ----------------
#>   Iter         Jac     Lambda      Eta     Dlt0     Dltn         Fnorm   Largest |f|
#>      0                                                    2.886812e+00  2.250000e+00
#>      1  N(9.6e-03) N            0.9544  10.1874   1.0187  5.787362e+05  1.070841e+03
#>      1             W   0.0932   0.9544   1.0187   1.0187  1.862204e+00  1.387696e+00
#>      2  B(1.9e-01) N            0.9658   0.6119   1.2239  1.003362e-01  4.468549e-01
#>      3  B(2.1e-01) N            0.6470   0.2521   0.5042  1.763418e-02  1.819978e-01
#>      4  B(1.3e-01) N            0.3486   0.2617   0.0598  5.952083e-02  3.020657e-01
#>      4             W   0.6248   0.3486   0.0598   0.1196  6.625459e-03  8.319100e-02
#>      5  B(1.5e-01) N            0.4212   0.1186   0.1186  3.368227e-03  5.868436e-02
#>      6  B(1.4e-01) N            0.9955   0.0177   0.0354  9.352918e-05  1.350724e-02
#>      7  B(1.8e-01) N            0.7560   0.0073   0.0146  9.330438e-06  3.933436e-03
#>      8  B(1.5e-01) N            0.8799   0.0020   0.0041  1.458619e-09  3.825434e-05
#>      9  B(1.5e-01) N            0.9969   0.0000   0.0000  1.640555e-13  4.501888e-07
#>     10  B(1.5e-01) N            0.9979   0.0000   0.0000  3.239033e-18  2.056389e-09
#> $x
#> [1] 1 1
#> 
#> $fvec
#> [1] 1.499777e-09 2.056389e-09
#> 
#> $termcd
#> [1] 1
#> 
#> $message
#> [1] "Function criterion near zero"
#> 
#> $scalex
#> [1] 1 1
#> 
#> $nfcnt
#> [1] 12
#> 
#> $njcnt
#> [1] 1
#> 
#> $iter
#> [1] 10
#> 

# final Broyden approximation of Jacobian (quite good)
z <- nleqslv(xstart, dslnex, jacobian=TRUE,control=list(btol=.01))
z$x
#> [1] 1 1
z$jac
#>           [,1]     [,2]
#> [1,] 1.9818090 2.005414
#> [2,] 0.9750673 3.007420
jacdsln(z$x)
#>      [,1] [,2]
#> [1,]    2    2
#> [2,]    1    3

# different initial start; not a very good final approximation
xstart <- c(0.5,2)
z <- nleqslv(xstart, dslnex, jacobian=TRUE,control=list(btol=.01))
z$x
#> [1] 1 1
z$jac
#>          [,1]     [,2]
#> [1,] 2.197044 2.182400
#> [2,] 3.210320 5.046055
jacdsln(z$x)
#>      [,1] [,2]
#> [1,]    2    2
#> [2,]    1    3

if (FALSE) { # \dontrun{
# no global strategy but limit stepsize
# but look carefully: a different solution is found
nleqslv(xstart, dslnex, method="Newton", global="none", control=list(trace=1,stepmax=5))

# but if the stepsize is limited even more the c(1,1) solution is found
nleqslv(xstart, dslnex, method="Newton", global="none", control=list(trace=1,stepmax=2))

# Broyden also finds the c(1,1) solution when the stepsize is limited
nleqslv(xstart, dslnex, jacdsln, method="Broyden", global="none", control=list(trace=1,stepmax=2))
} # }

# example with a singular jacobian in the initial guess
f <- function(x) {
    y <- numeric(3)
    y[1] <- x[1] + x[2] - x[1]*x[2] - 2
    y[2] <- x[1] + x[3] - x[1]*x[3] - 3
    y[3] <- x[2] + x[3] - 4
    return(y)
}

Jac <- function(x) {
    J <- matrix(0,nrow=3,ncol=3)
    J[,1] <- c(1-x[2],1-x[3],0)
    J[,2] <- c(1-x[1],0,1)
    J[,3] <- c(0,1-x[1],1)
    J
}

# exact solution
xsol <- c(-.5, 5/3 , 7/3)
xsol
#> [1] -0.500000  1.666667  2.333333

xstart <- c(1,2,3)
J <- Jac(xstart)
J
#>      [,1] [,2] [,3]
#> [1,]   -1    0    0
#> [2,]   -2    0    0
#> [3,]    0    1    1
rcond(J)
#> [1] 0

z <- nleqslv(xstart,f,Jac, method="Newton",control=list(trace=1,allowSingular=TRUE))
#>   Algorithm parameters
#>   --------------------
#>   Method: Newton  Global strategy: double dogleg (initial trust region = -2)
#>   Maximum stepsize = 1.79769e+308
#>   Scaling: fixed
#>   ftol = 1e-08 xtol = 1e-08 btol = 0.001 cndtol = 1e-12
#> 
#>   Iteration report
#>   ----------------
#>   Iter         Jac     Lambda      Eta     Dlt0     Dltn         Fnorm   Largest |f|
#>      0                                                    3.000000e+00  2.000000e+00
#>      1  N(0.0e+00) N            0.9535   1.2247   2.4495  2.500000e-01  5.000000e-01
#>      2  N(2.0e-01) N            0.8000   0.6124   1.2247  1.562500e-02  1.250000e-01
#>      3  N(2.6e-01) N            0.9660   0.1179   0.2357  2.999644e-28  1.731948e-14
all.equal(z$x,xsol)
#> [1] TRUE
```
