
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- genReadsForGene(gene, c(2/9, 3/9, 4/9))

matchIso(gene, reads=reads)

######################

reads2 <- genReadsForGene(gene, c(2/9, 3/9, 4/9), paired=TRUE,
                          normalInsert=c(100, 10))

m2 <- matchIso(gene, reads=reads2)
m2

reads3 <- reads2 
attr(reads3, "paired") <- FALSE

m3 <- matchIso(gene, reads=reads3)
m3

m <- cbind(matrix(m2, ncol=2, byrow=TRUE), 
           matrix(m3, ncol=2, byrow=TRUE))

all(apply(m, 1, function(x) {
  r <- lapply(strsplit(x, ""), as.numeric)
  all(r[[1]] == r[[2]]) && all(r[[1]] == (r[[3]] & r[[4]]) + 0)
}))
