
library(splicing)

options(width=60)

fileDir <- system.file("test_files", package="splicing")

checkorder <- function(x) {
  x <- matrix(x, nc=2, byrow=TRUE)
  ! is.unsorted(pmin(x[,1], x[,2]))
}

f1 <- readSAM(sprintf("%s/%s", fileDir, "single.sam"))
f2 <- readSAM(sprintf("%s/%s", fileDir, "paired.sam"))
f2$position
checkorder(f2$position)
f3 <- readSAM(sprintf("%s/%s", fileDir, "single.bam"))
f4 <- readSAM(sprintf("%s/%s", fileDir, "paired.bam"))
all(f4$position == f2$position)
f5 <- readSAM(sprintf("%s/%s", fileDir, "single-sorted.bam"))
f6 <- readSAM(sprintf("%s/%s", fileDir, "paired-sorted.bam"))
all(f6$position == f2$position)

r1 <- readSAM(sprintf("%s/%s", fileDir, "single-sorted.bam"),
              "seq1:200-400")
range(r1$position)

r2 <- readSAM(sprintf("%s/%s", fileDir, "paired-sorted.bam"),
              "seq1:200-400")
range(r2$position)
r2$pairpos

