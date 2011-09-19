
library(splicing)

options(width=60)

fileDir <- system.file("test_files", package="splicing")

r2 <- readSAM(sprintf("%s/%s", fileDir, "paired-sorted.bam"),
              "seq1:200-400")

f1 <- filterReads(r2, "noPaired")
f2 <- filterReads(r2, FILTER("noPaired"))
f3 <- filterReads(r2, "noPaired", FILTER("noPaired"))

r1 <- readSAM(sprintf("%s/%s", fileDir, "paired-sorted.bam"))
r3 <- filterReads(r1, FILTER("keepRegion", sequence="seq1",
                             start=200, end=400))

mapply(r2, r3, FUN=function(x,y) length(x)==length(y) && all(x==y))
all(sort(r2$qname) == sort(r3$qname))
all(sort(r2$attributes) == sort(r3$attributes))

