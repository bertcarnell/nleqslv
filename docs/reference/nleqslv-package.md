# The nleqslv package provides two algorithms for solving (dense) nonlinear systems of equations.

The methods provided are

- a Broyden Secant method where the matrix of derivatives is updated
  after each major iteration using the Broyden rank 1 update.

- a full Newton method where the Jacobian matrix of derivatives is
  recalculated at each iteration

Both methods utilize global strategies such as line search or trust
region methods whenever the standard Newton/Broyden step does not lead
to a point closer to a root of the equation system. Both methods can
also be used without a norm reducing global strategy. Line search may be
either cubic, quadratic or geometric. The trust region methods are
either the double dogleg method, the Powell single dogleg method or a
Levenberg-Marquardt type method.

## Details

There is a facility for specifying that the Jacobian is banded; this can
significantly speedup the calculation of a numerical Jacobian when the
number of sub- and super diagonals is small compared to the size of the
system of equations. For example the Jacobian of a tridiagonal system
can be calculated with only three evaluations of the function.

The package provides an option to attempt to solve the system of
equations when the Jacobian is singular or ill-conditioned using an
approximation to the Moore-Penrose pseudoinverse of the Jacobian.

The algorithms provided in this package are derived from Dennis and
Schnabel (1996). The code is written in Fortran 77 and Fortran 95 and
uses Lapack and BLAS routines as provided by the R system.

## References

Dennis, J.E. Jr and Schnabel, R.B. (1996), *Numerical Methods for
Unconstrained Optimization and Nonlinear Equations*, Siam.

## See also

Useful links:

- <https://bertcarnell.github.io/nleqslv/>

- <https://github.com/bertcarnell/nleqslv/>

- Report bugs at <https://github.com/bertcarnell/nleqslv/issues/>

## Author

**Maintainer**: Rob Carnell <bertcarnell@gmail.com>
([ORCID](https://orcid.org/0009-0009-0465-7564))

Authors:

- Berend Hasselman <bhh@xs4all.nl> (Original Author)
