
library(splicing)

check <- function(mean, var, numdevs, minlength) {
  
  start <- max(minlength, mean-numdevs*sqrt(var))
  end <- mean+numdevs*sqrt(var)
  f1 <- dnorm(start:end, mean, sqrt(var))
  f2 <- .Call("R_splicing_normal_fragment", as.double(mean), as.double(var),
              as.double(numdevs), as.integer(minlength),
              PACKAGE="splicing")
  (start == f2$fragmentStart && length(f1) == length(f2$fragmentProb) &&
   all(abs(f1 - f2$fragmentProb) < 1e-14))
}

check(100, 100, 4, 50)
check(250, 400, 4, 66)
check(50, 400, 4, 5)
