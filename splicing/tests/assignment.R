
library(splicing)

gene <- createGene(list(c(1,20), c(31,40)),
                   list(c(1), c(1,2)))

assignmentMatrix(gene, readLength=5, paired=TRUE,
                 fragmentStart=13, fragmentProb=(1:5)/sum(1:5),
                 normalMean=0, normalVar=0, numDevs=0)

