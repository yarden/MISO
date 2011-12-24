
library(splicing)

set.seed(42)
options(width=60)

gene <- createGene(list(c(1,100), c(201,300), c(401,500)),
                   list(c(1,2), c(1,3), c(1,2,3)))

ex <- constitutiveExons(gene, 10)
ex
ex$start
ex$end

gene <- createGene(list(c(1,500), c(101,500)), list(1,2))

ex <- constitutiveExons(gene, 10)
ex
ex$start
ex$end

ex <- constitutiveExons(gene, 10, mode="all")
ex
ex$start
ex$end

