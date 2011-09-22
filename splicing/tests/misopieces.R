
library(splicing)

options(width=60)

## MVRNORM

set.seed(42)
mu <- rep(1, 5)
sigma <- sqrt(2)
len <- as.integer(length(mu))
.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing")
.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing")
.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing")
.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing")
.Call("R_splicing_mvrnorm", mu, sigma, len, PACKAGE="splicing")

## LOGIT_INV

set.seed(42)
x <- runif(10)
.Call("R_splicing_logit_inv", x, as.integer(length(x)-1L), PACKAGE="splicing")

## SCORE_JOINT

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

assignment <- as.integer(sample(0:2, noReads(reads), replace=TRUE))
noreads <- as.integer(noReads(reads))
psi <- c(2/10, 3/10, 5/10)
hyper <- rep(1, noIso(gene))
effisolen <- pmax(isoLength(gene)[[1]] - getReadLength(reads) + 1, 0)
isoscores <- -log(effisolen)

.Call("R_splicing_score_joint", assignment, noreads, psi, hyper,
      as.integer(effisolen), isoscores, PACKAGE="splicing")

assignment <- as.integer(sample(0:2, noReads(reads), replace=TRUE,
                                prob=reads$sampleProb))

.Call("R_splicing_score_joint", assignment, noreads, psi, hyper,
      as.integer(effisolen), isoscores, PACKAGE="splicing")

## REASSIGN_SAMPLES

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

matches <- matchIso(gene, reads=reads)
match.order <- order(apply(matches, 2, paste, collapse=""))
psi <- c(2/10, 3/10, 5/10)
noiso <- 3L

.Call("R_splicing_reassign_samples", matches, match.order, psi, noiso,
      PACKAGE="splicing")

## DRIFT_PROPOSAL

set.seed(42)

## init
.Call("R_splicing_drift_proposal", 0L, NULL, NULL, 0.0, NULL, NULL, 2L,
      PACKAGE="splicing")
.Call("R_splicing_drift_proposal", 0L, NULL, NULL, 0.0, NULL, NULL, 3L,
      PACKAGE="splicing")

## propose
psi <- c(3/4, 1/4)
alpha <- 0
sigma <- 0.05
noiso <- length(psi)
.Call("R_splicing_drift_proposal", 1L, psi, alpha, sigma, NULL, NULL, noiso,
      PACKAGE="splicing")

psi <- c(1/3, 1/3, 1/3)
alpha <- c(1/2, 1/2)
sigma <- 0.05
noiso <- length(psi)
.Call("R_splicing_drift_proposal", 1L, psi, alpha, sigma, NULL, NULL, noiso,
      PACKAGE="splicing")

## score
psi <- c(3/4, 1/4)
alpha <- 0
sigma <- 0.05
otherpsi <- c(51/100, 49/100)
otheralpha <- 3/100
noiso <- length(psi)
.Call("R_splicing_drift_proposal", 2L, psi, alpha, sigma, otherpsi,
      otheralpha, noiso,
      PACKAGE="splicing")

## METROPOLIS_HASTINGS_RATIO

set.seed(42)
gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

matches <- matchIso(gene, reads=reads)
match.order <- order(apply(matches, 2, paste, collapse=""))-1L
psi <- c(2/10, 3/10, 5/10)
alpha <- c(1/3, 2/3)
noiso <- 3L

set.seed(42)
assignment <- .Call("R_splicing_reassign_samples", matches, match.order,
                    psi, noiso, PACKAGE="splicing")

noreads <- noReads(reads)
psinew <- c(5/10, 3/10, 2/10)
alphanew <- c(1/2, 1/2)
sigma <- 0.05
effisolen <- pmax(isoLength(gene)[[1]] - getReadLength(reads) + 1L, 0L)
hyperp <- rep(1, noiso)
isoscores <- -log(effisolen)

.Call("R_splicing_metropolis_hastings_ratio", assignment, noreads,
      psinew, alphanew, psi, alpha, sigma, noiso, as.integer(effisolen),
      hyperp, isoscores, 1L,
      PACKAGE="splicing")

.Call("R_splicing_metropolis_hastings_ratio", assignment, noreads,
      psi, alpha, psinew, alphanew, sigma, noiso, as.integer(effisolen),
      hyperp, isoscores, 1L,
      PACKAGE="splicing")

