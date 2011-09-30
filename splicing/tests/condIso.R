
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

geneComplexity(gene, readLength=33L)
geneComplexity(gene, readLength=33L, paired=TRUE,
               normalMean=166, normalVar=100, numDevs=4)

M1 <- assignmentMatrix(gene, readLength=33L)
M2 <- assignmentMatrix(gene, readLength=33L, paired=TRUE,
                       normalMean=166, normalVar=100, numDevs=4)
kappa(M1, exact=TRUE)
kappa(M2, exact=TRUE)

