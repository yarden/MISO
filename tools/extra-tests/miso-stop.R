
library(splicing)

mus <- readGFF3("/Users/gaborcsardi/tmp/Mus_musculus.gff3")

gene <- selectGenes(mus, 5)

percentAccepted <- function(x) {
  runData(x)$noAccepted / (runData(x)$noAccepted + runData(x)$noRejected)
}
norm <- function(x) x/sum(x)
exp <- norm((10:1)^4)
reads <- simulateReads(gene, readLength=35, noReads=4000, expression=exp)
rexp <- realPsi(gene, reads=reads)

mres <- MISO(gene, reads=reads, noIterations=5000, noBurnIn=2500, noLag=10,
             stopCond="convMean", noChains=2)
mresln <- MISO(gene, reads=reads, noIterations=5000, noBurnIn=2500,
               noLag=10, start="linear", noChains=2, stopCond="convMean")

cor(exp, postMean(mres))
cor(exp, postMean(mresln))
sum(abs(exp-postMean(mres)))
sum(abs(exp-postMean(mresln)))
percentAccepted(mres)
percentAccepted(mresln)

####################

gene3 <- selectGenes(mus, 3)

exp3 <- norm((3:1)^2)
reads3 <- simulateReads(gene3, readLength=35, noReads=4000, expression=exp3)
rexp3 <- realPsi(gene3, reads=reads3)

mres3 <- MISO(gene3, reads=reads3, noIterations=5000, noChains=2,
              noBurnIn=500, noLag=10, stopCond="convMean")
mresln3 <- MISO(gene3, reads=reads3, noIterations=5000, noChains=2,
              noBurnIn=500, noLag=10, start="linear", stopCond="convMean")

cor(exp3, postMean(mres3))
cor(exp3, postMean(mresln3))
sum(abs(exp3-postMean(mres3)))
sum(abs(exp3-postMean(mresln3)))
percentAccepted(mres3)
percentAccepted(mresln3)

