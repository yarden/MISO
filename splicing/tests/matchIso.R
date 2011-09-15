
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- simulateReads(gene, expression=c(2/9, 3/9, 4/9), noReads=100L,
                       readLength=33L)

matchIso(gene, reads=reads)

######################

reads2 <- simulateReads(gene, expression=c(2/10, 3/10, 5/10), noReads=100L,
                        paired=TRUE, readLength=33L, normalMean=166,
                        normalVar=100, numDevs=4)

m2 <- matchIso(gene, reads=reads2, paired=FALSE)

m3 <- matchIso(gene, reads=reads2,
               normalMean=166, normalVar=100, numDevs=4)
m3

all(sapply(1:ncol(m3[[1]]),
           function(x) m3[[1]][,x][ getIsoform(reads2)[2*x-1]+1 ]) != 0)
all(sapply(1:ncol(m3[[1]]),
           function(x) m3[[1]][,x][ getIsoform(reads2)[2*x]+1 ]) != 0)

mean(m3[[2]] [ m3[[1]] != 0 ])

M <- cbind(apply(m2[,seq(1,ncol(m2),by=2)], 2, paste, collapse=""),
           apply(m2[,seq(2,ncol(m2),by=2)], 2, paste, collapse=""),
           apply((m3[[1]] != 0) + 0, 2, paste, collapse=""))

# Any contradiction?

which(apply(M, 1, function(a) {
  b <- lapply(strsplit(a, split=""), as.numeric)
  any((!b[[1]] & b[[3]]) | (!b[[2]] & b[[3]]))
}))
