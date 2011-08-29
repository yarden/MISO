
library(splicing)

fileDir <- system.file("test_files", package="splicing")

f1 <- readSAM(sprintf("%s/%s", fileDir, "single.sam"))
f2 <- readSAM(sprintf("%s/%s", fileDir, "paired.sam"))
f3 <- readSAM(sprintf("%s/%s", fileDir, "single.bam"))
f4 <- readSAM(sprintf("%s/%s", fileDir, "paired.bam"))
f5 <- readSAM(sprintf("%s/%s", fileDir, "single-sorted.bam"))
f6 <- readSAM(sprintf("%s/%s", fileDir, "paired-sorted.bam"))

r1 <- readSAM(sprintf("%s/%s", fileDir, "single-sorted.bam"),
              "seq1:200-400")
range(r1$position)
r2 <- readSAM(sprintf("%s/%s", fileDir, "paired-sorted.bam"),
              "seq1:200-400")
range(r2$position)
