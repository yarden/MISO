
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

set.seed(42)
mres1 <- MISO(gene, reads=reads, readLength=35L)

mm <- matchIso(gene, reads=reads)
il <- isoLength(gene)[[1]]

set.seed(42)
mres2 <- MISO.Trinity(mm, isoLength=il, readLength=35L)

all(postMean(mres1) == postMean(mres2))

###############

set.seed(42)
reads2 <- simulateReads(gene, expression=c(2/10, 3/10, 5/10), paired=TRUE,
                        noReads=1000L, readLength=33, normalMean=116,
                        normalVar=50, numDevs=4)

set.seed(42)
dres1 <- MISO(gene, reads=reads2, readLength=33L, normalMean=116, paired=TRUE,
              fragmentStart=0L, normalVar=50, numDevs=4)

mm <- matchIso(gene, reads=reads2, paired=TRUE, normalMean=116, normalVar=50,
               numDevs=4)
il <- isoLength(gene)[[1]]

set.seed(42)
dres2 <- MISO.Trinity(mm[[1]], fragmentLength=mm[[2]], il, readLength=33L,
                      normalMean=116, paired=TRUE, fragmentStart=0L,
                      normalVar=50, numDevs=4)

all(postMean(dres1) == postMean(dres2))

