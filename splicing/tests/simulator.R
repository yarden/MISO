
library(splicing)

set.seed(42)
options(width=60)

## Single end, TODO

gene <- createGene(list(c(1,1000), c(2001,3000)),
                   list(c(1), c(1,2)))

reads <- simulateReads(gene, 1, expression=c(1,0), noReads=5000,
                       readLength=250, paired=FALSE)

## Paired end

gene <- createGene(list(c(1,20), c(31,40)), list(c(1), c(1,2)))

reads <- simulateReads(gene, 1, expression=c(1,0), noReads=100,
                       readLength=5, paired=TRUE,
                       fragmentProb=(1:5)/sum(1:5), fragmentStart=13,
                       normalMean=0, normalVar=0, numDevs=0)

reads$sampleProb

reads <- simulateReads(gene, 1, expression=c(1,0), noReads=100,
                       readLength=5, paired=TRUE,
                       fragmentProb=(1:5)/sum(1:5), fragmentStart=18,
                       normalMean=0, normalVar=0, numDevs=0)

reads$sampleProb

## 

gene <- createGene(list(c(1,1000), c(2001,3000)),
                   list(c(1), c(1,2)))

reads <- simulateReads(gene, 1, expression=c(1,0), noReads=5000,
                       readLength=30, paired=TRUE, normalMean=250,
                       normalVar=400, numDevs=4)

hist(reads$position[seq(1,10000,by=2)], 16, plot=FALSE)$count

reads$paired && all(reads$isoform == 0) && all(reads$cigar=="30M")

