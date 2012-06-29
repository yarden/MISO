
library(splicing)

set.seed(42)
options(width=60)

##################################################
## A gene for which the bias cannot be estimated

gene <- createGene(list(c(1,100), c(1,300), c(201,300)),
                   list(c(1,3), 2))

assignmentMatrix(gene, gene=1L, readLength=30)
cbind("01"=c(0,200-71),
      "10"=c(100-71,0),
      "11"=c(71-0+171-100,71-0+271-200))

assignmentMatrix(gene, gene=1L, readLength=30, bias=1L)
cbind("01"=c(0,200^2-71^2),
      "10"=c(100^2-71^2,0),
      "11"=c(71^2-0^2+171^2-100^2,71^2-0^2+271^2-200^2))

assignmentMatrix(gene, gene=1L, readLength=30, bias=2L)
cbind("01"=c(0,200^3-71^3),
      "10"=c(100^3-71^3,0),
      "11"=c(71^3-0^3+171^3-100^3,71^3-0^3+271^3-200^3))

reads <- simulateReads(gene, readLength=30, noReads=2000,
                       expression=c(0.7,0.2,0.1))

s1 <- solveIso(gene, reads=reads)
tryres <- try(s2 <- solveIsoLinBias(gene, reads=reads), silent=TRUE)
inherits(tryres, "try-error")

#################################################
## Linear bias can be estimated, but not the
## quadratic bias

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(2,3), 1, 2, c(1,2,3)))

reads <- simulateReads(gene, readLength=30, noReads=2000,
                       expression=c(0.2, 0.3, 0.3, 0.1, 0.1, 0))

s1 <- solveIso(gene, reads=reads)
s2 <- solveIsoLinBias(gene, reads=reads)
tryres <- try(s3 <- solveIsoQuadBias(gene, reads=reads), silent=TRUE)
inherits(tryres, "try-error")

round(s1$expression, 4)
round(s2$expression, 4)
round(s2$a, 4)

#################################################
## Quadratic bias can be estimated here

gene <- createGene(list(c(1,100), c(201,300), c(401,500), c(601,700)),
                   list(c(1,2), c(1,3), c(2,3), c(2,4)))

reads <- simulateReads(gene, readLength=30, noReads=200,
                       expression=c(0.2, 0.3, 0.3, 0.1))

s1 <- solveIso(gene, reads=reads)
s2 <- solveIsoLinBias(gene, reads=reads)
s3 <- solveIsoQuadBias(gene, reads=reads)

round(s1$expression, 4)
round(s2$expression, 4)
round(s3$expression, 4)
round(s2$a, 4)
round(s3$a, 4)
round(s3$b, 4)
