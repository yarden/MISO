
## Generate, match and estimate overlapping paired-end reads

library(splicing)

dist <- .Call("R_splicing_normal_fragment", 50, 100, 4, 35L,
              PACKAGE="splicing")
dist$fragmentStart

## 

gene <- createGene(list(c(1,1000), c(2001, 3000)),
                   list(c(1,2), c(1)))

rl <- 50
normalMean <- 75
normalVar <-  75

reads <- simulateReads(gene, expression=c(2/10, 8/10),
                       noReads=1000L, readLength=rl, paired=TRUE,
                       normalMean=normalMean, normalVar=normalVar,
                       numDevs=4)

## 

matches <- matchIso(gene, reads=reads, normalMean=normalMean,
                    normalVar=normalVar, numDevs=4)

any(colSums(matches[[1]] !=  0) == 0)
any(colSums(matches[[2]] != -1) == 0)

##

misores <- MISO(gene, reads=reads, normalMean=normalMean,
                normalVar=normalVar, numDevs=4)
any(colSums(misores$matchMatrix != 0) == 0)

##

solveres <- solveIso(gene, reads=reads, normalMean=normalMean,
                     normalVar=normalVar, numDevs=4)
any(colSums(solveres$match != 0) == 0)

