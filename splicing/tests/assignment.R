
library(splicing)

gene <- createGene(list(c(1,20), c(31,40)),
                   list(c(1), c(1,2)))

assignmentMatrix(gene, readLength=5L, paired=TRUE,
                 fragmentStart=20L, fragmentProb=1,
                 normalMean=15, normalVar=0, numDevs=0)

###

gene <- createGene(list(c(1,10), c(21,30), c(41,50)),
                   list(c(1,2), c(1,3), c(1,2,3)))

assignmentMatrix(gene, readLength=5L, paired=TRUE,
                 fragmentStart=15L, fragmentProb=1,
                 normalMean=15, normalVar=0, numDevs=0)

###

gene <- createGene(list(c(1,20), c(31,40)),
                   list(c(1), c(1,2)))

assignmentMatrix(gene, readLength=5, paired=TRUE,
                 fragmentStart=13, fragmentProb=(1:5)/sum(1:5),
                 normalMean=0, normalVar=0, numDevs=0)


###

