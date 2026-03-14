# Copyright (c) Rob Carnell 2026

# More, Garbow, Hillstrom: Testing Unconstrained Optimization Software
# ACM Trans. Math. Software, 7, March 1981, 17--41
# Function as shown in this paper is for optimizing not for non linear equations
# from the paper it is not clear what was used by the authors for their tests
# of non linear equation solvers

# function 27
# Brown almost linear function

brown <- function(x) {
  n <- length(x)
  y <- numeric(n)

  y[1:(n-1)] <- x[1:(n-1)] + sum(x[1:n]) - (n + 1)
  y[n] <- prod(x[1:n]) - 1.0

  y
}

brown_xstart <- function(n) rep(1,n)/2

# function 28
# Discrete boundary value problem

dcbval <- function (x) {
  n <- length(x)
  y <- numeric(n)

  h <- 1/(n+1)

  y <- 2*x + .5*h^2 * (x+1+h*c(1:n))^3
  y[2:n]     <- y[2:n]     - x[1:(n-1)]
  y[1:(n-1)] <- y[1:(n-1)] - x[2:n]

  y
}

dcbval_xstart <- function(n) c(1:n) * (c(1:n)-n-1)/(n+1)^2

# function 29
# Discrete integral equation

dciequ <- function(x) {
  n <- length(x)
  y <- numeric(n)

  h <- 1/(n+1)

  # do it the stupid (and slow?) way

  for (k in 1:n) {
    tk <-  k/(n+1)

    sum1 <- 0
    for (j in 1:k) {
      tj <- j*h
      sum1 <- sum1 + tj * (x[j] + tj + 1.0)^3
    }

    # do this with while since R will count down with for loop
    sum2 <- 0
    j <- k+1
    while(j <= n) {
      tj  <- j*h
      sum2 <- sum2 + (1.0 - tj) * (x[j] + tj + 1.0)^3
      j <- j + 1
    }
    y[k] = x[k] + h*((1.0 - tk) * sum1 + tk * sum2 ) / 2.0
  }

  y
}

dciequ_xstart <- function(n) c(1:n) * (c(1:n)-1-n)/(n+1)^2

# function 7
# Helical valley function

helval <- function(x) {
  c1 <- 0.25
  c2 <- 0.5

  stopifnot(length(x) == 3)
  y <- numeric(3)

  tpi <- 8*atan(1)
  if (x[1] > 0)
    temp1 <- atan(x[2]/x[1])/tpi
  else if (x[1] < 0)
    temp1 <- atan(x[2]/x[1])/tpi + c2
  else
    temp1 <- sign(c1,x[2])

  y[1] <- 10*(x[3] - 10*temp1)
  y[2] <- 10*(sqrt(x[1]^2+x[2]^2) - 1)
  y[3] <- x[3]

  y
}

helval_xstart <- c(-1,0,0)

# Powell singular function
# function 13

pwlsng <- function(x) {
  stopifnot(length(x) == 4)
  y <- numeric(4)
  y[1] <- x[1] + 10*x[2]
  y[2] <- sqrt(5)*(x[3] - x[4])
  y[3] <- (x[2] - 2*x[3])^2
  y[4] <- sqrt(10)*(x[1] - x[4])^2

  y
}

pwlsng_xstart <- c(3, -1, 0, 1)

# Rosenbrock function
# function 1

rosbrk <- function(x) {
  stopifnot(length(x) == 2)
  y <- numeric(2)

  y[1] <- 1 - x[1]
  y[2] <- 10*(x[2] - x[1]^2)

  y
}

rosbrk_xstart <- c(-1.2, 1)

# vardim function
# Function 25
# The function here comes from
#  http://www.zib.de/Numerik/numsoft/CodeLib/codes/newton_testset/problems/VariablyDim.f

# known solutions
vardimsol <- function(n) {
  x <- numeric(n)
  x[1:n] <- 1
  x
}

vardim <- function(x) {
  n <- length(x)
  s <- sum(c(1:n)*(x-1))
  y <- x - 1 + s * (1+2*s^2) * c(1:n)
}

vardiminit <- function (n) {
  x <- (1 - c(1:n) / n)
}

# Watson function
# code adapted from test problems f90 for Minpack
# very difficult function
# function 20

watson <- function(x) {

  n <- length(x)
  y <- numeric(n)

  for(i in 1:29) {
    ti <- i/29

    sum1 <- 0
    t    <- 1
    for (j in 2:n) {
      sum1 <- sum1 + (j-1) * t * x[j]
      t    <- ti * t
    }

    sum2 <- 0
    t    <- 1
    for (j in 1:n) {
      sum2 <- sum2 + t * x[j]
      t    <- ti * t
    }

    t <- 1 / ti
    for(k in 1:n) {
      y[k] <- y[k] + t * (sum1-sum2^2-1) * (k-1-2*ti*sum2)
      t <- t * ti
    }
  }

  y[1] = y[1] + 3*x[1] - 2*x[1]*x[2] + 2*x[1]^3
  y[2] = y[2] + x[2] - x[1]^2 - 1

  y
}

watson_xstart <- function(n) numeric(n)

# Wood function
# function 14

wood <- function(x) {
  c1 <- 2.0e2
  c2 <- 2.02e1
  c3 <- 1.98e1
  c4 <- 1.8e2

  stopifnot(length(x) == 4)
  y <- numeric(4)

  y[1] <- -c1*x[1]*(x[2] - x[1]^2) - (1 - x[1])
  y[2] <- c1*(x[2] - x[1]^2) + c2*(x[2] - 1) + c3*(x[4]-1)
  y[3] <- -c4*x[3]*(x[4] - x[3]^2) - (1 - x[3])
  y[4] <- c4*(x[4] - x[3]^2) + c2*(x[4] - 1) + c3*(x[2]-1)

  y
}

wood_xstart <- c(-3,-1,-3,-1)

# Powell Badly scaled function
# function 3

pwlbsc <- function(x) {
  c1 <- 1.e4
  c2 <- 1.0001

  stopifnot(length(x) == 2)
  y <- numeric(2)
  y[1] <- c1*x[1]*x[2] - 1
  y[2] <- exp(-x[1]) + exp(-x[2]) - c2

  y
}

pwlbscjac <- function(x) {
  stopifnot(length(x) == 2)
  c1 <- 1.e4
  c2 <- 1.0001
  J <- matrix(NA, 2, 2)
  J[1,1] <- c1*x[2]
  J[1,2] <- c1*x[1]
  J[2,1] <- -exp(-x[1])
  J[2,2] <- -exp(-x[2])
  J
}

pwlbsc_xstart <- c(0,1)

# Dennis & Schnabel,1996,"Numerical methods for unconstrained optimization and nonlinear equations", SIAM
# example 6.5.1 page 149

dslnex <- function(x) {
  stopifnot(length(x) == 2)
  y <- numeric(2)
  y[1] <- x[1]^2 + x[2]^2 - 2
  y[2] <- exp(x[1]-1) + x[2]^3 - 2
  y
}

jacdsln <- function(x) {
  stopifnot(length(x) == 2)
  Df <- matrix(NA, 2, 2)
  Df[1,1] <- 2*x[1]
  Df[1,2] <- 2*x[2]
  Df[2,1] <- exp(x[1]-1)
  Df[2,2] <- 3*x[2]^2

  Df
}

# give parameters names which must appear in results
dslnex_xstart <- c(x1=1.5, x2=2)

# Nocedal & Wright (2006) example 11.31

f_nocedal <- function(x) {
  stopifnot(length(x) == 2)
  y <- numeric(2)
  y[1] <- (x[1]+3)*(x[2]^2-7) + 18
  y[2] <- sin(x[2]*exp(x[1])-1)
  y
}

jac_nocedal <- function(x) {
  stopifnot(length(x) == 2)
  J <- matrix(0, nrow=2, ncol=2)
  J[1,1] <- x[2]^2 - 7
  J[1,2] <- (x[1] + 3)*(2*x[2])
  J[2,2] <- cos(x[2]*exp(x[1])-1)*exp(x[1])
  J[2,1] <- x[2]*J[2,2]  # cos(x[2]*exp(x[1])-1)*x[2]*exp(x[1])
  J
}

nocedal_xstart <- c(-1, 2)
