
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/9, 3/9, 4/9), noReads=100L,
                       readLength=33L)

matchIso(gene, reads=reads)

######################

reads2 <- simulateReads(gene, expression=c(2/10, 3/10, 4/10), noReads=100L,
                        paired=TRUE, readLength=33L, normalMean=166,
                        normalVar=100, numDevs=4)

m3 <- matchIso(gene, reads=reads2, normalMean=166, normalVar=100, numDevs=4)
m3
