
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

reads <- genReadsForGene(gene, c(2/9, 3/9, 4/9))

condIso(gene, mode="matrix")
