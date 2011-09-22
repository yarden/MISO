
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

mres <- MISO(gene, reads=reads, readLength=35L)

postMean(mres)

reads2 <- simulateReads(gene, expression=c(2/10, 3/10, 5/10), paired=TRUE,
                        noReads=1000L, readLength=33, normalMean=116,
                        normalVar=50, numDevs=4)

mres2 <- MISO(gene, reads=reads2, readLength=33L, normalMean=116, paired=TRUE,
              fragmentStart=0L, normalVar=50, numDevs=4)

postMean(mres2)

## Multiple chains

mres3 <- MISO(gene, reads=reads, noChains=9, noIterations=1000)

postMean(mres3)

## Multiple chains, different (random) start positions

mres4 <- MISO(gene, reads=reads, noChains=9, noIterations=1, noBurnIn=0,
              noLag=1, start="random")
mres4$samples

