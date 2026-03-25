# Copyright (c) 2026 Rob Carnell

################################################################################
#
# Create R/nleqslv-iterationreport.R
#
################################################################################

require(devtools)
devtools::load_all()

output_file <- file.path("R", "nleqslv-iterationreport.R")
input_file <- file.path("etc", "nleqslv-iterationreport-template.txt")

X <- readLines(con = input_file)

# Dennis & Schnabel,1996,"Numerical methods for unconstrained optimization and nonlinear equations", SIAM
# example 6.5.1 page 149

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

xstart <- c(2,0.5)

table_linesearch <- capture.output({
  nleqslv(xstart, dslnex, global="qline", control=list(trace=1))
})[10:29]

table_dbldog <- capture.output({
  nleqslv(xstart, dslnex, global="dbldog", jacobian=TRUE,
          control=list(trace=1, delta="cauchy"))
})[10:24]

table_pwldog <- capture.output({
  nleqslv(xstart, dslnex, global="pwldog", jacobian=TRUE,
          control=list(trace=1, delta="cauchy"))
})[10:25]

table_hook <- capture.output({
  nleqslv(xstart, dslnex, global="hook", jacobian=TRUE,
          control=list(trace=1, delta="cauchy"))
})[10:21]

Xout <- c(X[1:(grep("<<table-linesearch>>", X) - 1)], table_linesearch)
Xout <- c(Xout, X[(grep("<<table-linesearch>>", X) + 1):(grep("<<table-dbldog>>", X) - 1)], table_dbldog)
Xout <- c(Xout, X[(grep("<<table-dbldog>>", X) + 1):(grep("<<table-pwldog>>", X) - 1)], table_pwldog)
Xout <- c(Xout, X[(grep("<<table-pwldog>>", X) + 1):(grep("<<table-hook>>", X) - 1)], table_hook)
Xout <- c(Xout, X[(grep("<<table-hook>>", X) + 1):length(X)])

Xout <- paste("#'", Xout)
Xout <- c("# DO NOT EDIT - use etc/nleqslv-iterationreport-creator.R",
          "",
          "# Copyright (c) 2026 Rob Carnell",
          "",
          Xout,
          "NULL")

writeLines(Xout, output_file)

