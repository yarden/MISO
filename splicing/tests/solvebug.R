
library(splicing)

gff <- system.file("test_files/ENSMUSG00000026173.gff", package="splicing")
readfile <- system.file("test_files/reads-4-10.Rdata", package="splicing")

gene <- readGFF3(gff)
load(readfile)
reads2 <- selectReads(reads, grep("read-5-", reads$qname))

res <- solveIso(gene, reads=reads2, paired=TRUE, normalMean=250+33+33,
                normalVar=250+33+33, numDevs=4)

mat <- matchIso(gene, reads=reads2, normalMean=250+33+33, normalVar=250+33+33,
                numDevs=4)

ass <- assignmentMatrix(gene, readLength=33,
                        paired=TRUE, normalMean=250+33+33,
                        normalVar=250+33+33, numDevs=4)

