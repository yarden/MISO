
library(splicing)
library(snow)

genelist <- lapply(seq(10,100,by=10), function(x) {
  createGene(exons=list(c(1,1000), c(2001, 2000+x), c(4001,5000)),
             isoforms=list(c(1,2,3), c(1,3)), id=sprintf("gene-%d", x))
})

genes <- do.call(merge, genelist)

set.seed(42)
readlist <- lapply(seq_along(genelist), function(x)
                   simulateReads(geneStructure=genes, gene=x, 
                                 qname=sprintf("read-%i-%%i", x),
                                 expression=c(2/10, 8/10),
                                 noReads=1000L, readLength=35))

reads <- do.call(mergeReads, readlist)

samfile <- paste(tempfile(), sep="", ".sam")
bamfile <- paste(tempfile(), sep="", ".bam")
bamfile2 <- tempfile()
writeSAM(reads, samfile)
SAMFile2BAMFile(samfile, bamfile)
sortBAMFile(bamfile, bamfile2)
indexBAMFile(paste(bamfile2, sep="", ".bam"))

cl <- makeCluster(10)

runMISO(genes, paste(bamfile2, sep="", ".bam"), snowCluster=cl, seed=42,
        verbose=FALSE)

stopCluster(cl)
