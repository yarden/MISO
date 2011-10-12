
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

## Check the pairs

load(system.file("test_files/mymm.Rdata", package="splicing"))

expr <- rep(1/noIso(mymm)[2], noIso(mymm)[2])

set.seed(42)
res <- .Call("R_splicing_simulate_paired_reads", mymm, 2L, expr, 173L, 33L, 
             NULL, 0L, 250+33+33, 250+33+33, 4, PACKAGE="splicing")

set.seed(42)
reads <- simulateReads(mymm, 2, expr, noReads=346/2, readLength=33,
                       paired=TRUE, normalMean=250+33+33, normalVar=250+33+33,
                       numDevs=4)

all(res$position == reads$position)
all(reads$position[reads$mypair+1] == reads$pairpos)
