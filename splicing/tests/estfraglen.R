
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,1000), c(2001,3000), c(4001,5000)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=rep(1,noIso(gene)),
                       noReads=1000, readLength=35, paired=TRUE,
                       normalMean=200, normalVar=200, numDevs=4)

tmp <- tempfile(fileext=c(".sam", ".bam", ""))
writeSAM(reads, tmp[1])
SAMFile2BAMFile(tmp[1], tmp[2])
sortBAMFile(tmp[2], tmp[3])
indexBAMFile(paste(tmp[3], sep="", ".bam"))

ex <- constitutiveExons(gene, 600)
es1 <- estimateFragLength(ex, readsfile=paste(tmp[3], sep="", ".bam"))
es1

# es2 <- estimateFragLength(ex, reads=reads)

