
library(splicing)
options(width=60)

load("~/tmp/mymm.Rdata")
geneStructure <- mymm

get.region <- function(gene) {
  start <- geneStructure$start[geneStructure$gid[gene]+1]
  end <- geneStructure$end[geneStructure$gid[gene]+1]
  seqid <- geneStructure$seqid[gene]
  seqid_str <- geneStructure$seqid_str[seqid+1]
  paste(sep="", seqid_str, ":", start, "-", end)
}

reads <- readSAM("/Users/gaborcsardi/tmp/reads-1-10_sorted.bam",
                 region=get.region(5))

mat <- matchIso(mymm, gene=5, reads=reads,
                paired=TRUE, normalMean=250+33+33, normalVar=250+33+33,
                numDevs=4)

ass <- assignmentMatrix(mymm, gene=5, paired=TRUE, readLength=33,
                        normalMean=250+33+33, normalVar=250+33+33,
                        numDevs=4)

ass

# sol1 <- solveIso(mymm, gene=5, reads=reads, paired=TRUE,
#                  normalMean=250+33+33, normalVar=250+33+33, numDevs=4)

## 

reads2 <- readSAM("/Users/gaborcsardi/tmp/reads-1-10_sorted.bam",
                 region=get.region(6))

mat2 <- matchIso(mymm, gene=6, reads=reads2,
                 paired=TRUE, normalMean=250+33+33, normalVar=250+33+33,
                 numDevs=4)

ass2 <- assignmentMatrix(mymm, gene=6, paired=TRUE, readLength=33,
                         normalMean=250+33+33, normalVar=250+33+33,
                         numDevs=4)

ass2

# sol2 <- solveIso(mymm, gene=6, reads=reads2,
#         paired=TRUE, normalMean=250+33+33, normalVar=250+33+33, numDevs=4)

