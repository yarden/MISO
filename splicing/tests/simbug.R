
library(splicing)

gff <- system.file("test_files/simbug.gff", package="splicing")
g <- readGFF3(gff)

set.seed(42)
tmp <- simulateReads(g, expression=c(0,0,1,0), readLength=33,
                     noReads=100, paired=TRUE, normalMean=316,
                     normalVar=316, numDevs=4)

grep("-", tmp$cigar, fixed=TRUE)

