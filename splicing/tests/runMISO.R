
library(splicing)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

set.seed(42)
reads <- simulateReads(gene, expression=c(2/10, 3/10, 5/10),
                       noReads=1000L, readLength=35)

samfile <- paste(tempfile(), sep="", ".sam")
bamfile <- paste(tempfile(), sep="", ".bam")
bamfile2 <- tempfile()
writeSAM(reads, samfile)
SAMFile2BAMFile(samfile, bamfile)
sortBAMFile(bamfile, bamfile2)
indexBAMFile(paste(bamfile2, sep="", ".bam"))

set.seed(42)
runMISO(gene, paste(bamfile2, sep="", ".bam"))

