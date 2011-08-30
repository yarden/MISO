
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=100L, readLength=20)

mres <- solveIso(gene, reads=reads, readLength=20L)
mres$expression

reads2 <- simulateReads(gene, expression=c(2/10, 3/10, 5/10), paired=TRUE,
                        noReads=100L, readLength=20, normalMean=90,
                        normalVar=50, numDevs=4)

mres2 <- solveIso(gene, reads=reads2, readLength=20L, paired=TRUE,
                  normalMean=90, normalVar=50, numDevs=4)
mres2$expression

