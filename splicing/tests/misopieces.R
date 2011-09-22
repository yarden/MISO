
library(splicing)

options(width=60)

## MVRNORM

set.seed(42)
mu <- cbind(rep(1, 5))
sigma <- sqrt(2)
len <- as.integer(length(mu))
as.vector(.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing"))
as.vector(.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing"))
as.vector(.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing"))
as.vector(.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing"))
as.vector(.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing"))

## LOGIT_INV

set.seed(42)
x <- cbind(runif(10))
as.vector(.Call("R_splicing_logit_inv", x, as.integer(length(x)-1L), 1L,
                PACKAGE="splicing"))

## SCORE_JOINT

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

assignment <- cbind(as.integer(sample(0:2, noReads(reads), replace=TRUE)))
noreads <- as.integer(noReads(reads))
psi <- cbind(c(2/10, 3/10, 5/10))
hyper <- rep(1, noIso(gene))
effisolen <- pmax(isoLength(gene)[[1]] - getReadLength(reads) + 1, 0)
isoscores <- -log(effisolen)

.Call("R_splicing_score_joint", assignment, noreads, 1L, psi, hyper,
      as.integer(effisolen), isoscores, PACKAGE="splicing")

assignment <- cbind(as.integer(sample(0:2, noReads(reads), replace=TRUE,
                                      prob=reads$sampleProb)))

.Call("R_splicing_score_joint", assignment, noreads, 1L, psi, hyper,
      as.integer(effisolen), isoscores, PACKAGE="splicing")

## REASSIGN_SAMPLES

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

matches <- matchIso(gene, reads=reads)
match.order <- order(apply(matches, 2, paste, collapse=""))-1L
psi <- cbind(c(2/10, 3/10, 5/10))
noiso <- 3L

as.vector(.Call("R_splicing_reassign_samples", matches, match.order, psi,
                noiso, 1L, PACKAGE="splicing"))

## DRIFT_PROPOSAL

set.seed(42)

## init
.Call("R_splicing_drift_proposal", 0L, NULL, NULL, 0.0, NULL, NULL, 2L, 1L,
      0L, NULL, NULL, PACKAGE="splicing")
.Call("R_splicing_drift_proposal", 0L, NULL, NULL, 0.0, NULL, NULL, 3L, 1L,
      0L, NULL, NULL, PACKAGE="splicing")

## propose
psi <- cbind(c(3/4, 1/4))
alpha <- cbind(0)
sigma <- 0.05
noiso <- nrow(psi)
.Call("R_splicing_drift_proposal", 1L, psi, alpha, sigma, NULL, NULL, noiso,
      1L, 0L, NULL, NULL, PACKAGE="splicing")

psi <- cbind(c(1/3, 1/3, 1/3))
alpha <- cbind(c(1/2, 1/2))
sigma <- 0.05
noiso <- nrow(psi)
.Call("R_splicing_drift_proposal", 1L, psi, alpha, sigma, NULL, NULL, noiso,
      1L, 0L, NULL, NULL, PACKAGE="splicing")

## score
psi <- cbind(c(3/4, 1/4))
alpha <- cbind(0)
sigma <- 0.05
otherpsi <- cbind(c(51/100, 49/100))
otheralpha <- cbind(3/100)
noiso <- length(psi)
.Call("R_splicing_drift_proposal", 2L, psi, alpha, sigma, otherpsi,
      otheralpha, noiso, 1L, 0L, NULL, NULL,
      PACKAGE="splicing")

## METROPOLIS_HASTINGS_RATIO

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

matches <- matchIso(gene, reads=reads)
match.order <- order(apply(matches, 2, paste, collapse=""))-1L
psi <- cbind(c(2/10, 3/10, 5/10))
alpha <- cbind(c(1/3, 2/3))
noiso <- 3L

set.seed(42)
assignment <- .Call("R_splicing_reassign_samples", matches, match.order,
                    psi, noiso, 1L, PACKAGE="splicing")

noreads <- noReads(reads)
psinew <- cbind(c(5/10, 3/10, 2/10))
alphanew <- cbind(c(1/2, 1/2))
sigma <- 0.05
effisolen <- pmax(isoLength(gene)[[1]] - getReadLength(reads) + 1L, 0L)
hyperp <- rep(1, noiso)
isoscores <- -log(effisolen)

.Call("R_splicing_metropolis_hastings_ratio", assignment, noreads, 1L,
      psinew, alphanew, psi, alpha, sigma, noiso, as.integer(effisolen),
      hyperp, isoscores, 1L,
      PACKAGE="splicing")

.Call("R_splicing_metropolis_hastings_ratio", assignment, noreads, 1L,
      psi, alpha, psinew, alphanew, sigma, noiso, as.integer(effisolen),
      hyperp, isoscores, 1L,
      PACKAGE="splicing")

## RNG_GET_DIRICHLET

set.seed(42)
alpha <- c(1,1,1,1)
alpha2 <- c(1,2,3,4)
.Call("R_splicing_rng_get_dirichlet", alpha, PACKAGE="splicing")
.Call("R_splicing_rng_get_dirichlet", alpha2, PACKAGE="splicing")
