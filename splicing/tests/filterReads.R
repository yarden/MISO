
library(splicing)

options(width=60)

fileDir <- system.file("test_files", package="splicing")

r2 <- readSAM(sprintf("%s/%s", fileDir, "paired-sorted.bam"),
              "seq1:200-400")

f1 <- filterReads(r2, "noPaired")
f2 <- filterReads(r2, FILTER("noPaired"))
f3 <- filterReads(r2, "noPaired", FILTER("noPaired"))

