
## MISO single chain vs multiple chains, for a gene with many
## isoforms and high complexity

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

mres <- MISO(gene, reads=reads, noIterations=5000, noBurnIn=2500, noLag=10)
mres2 <- MISO(gene, reads=reads, noChains=5, noIterations=1000,
              noBurnIn=500, noLag=10)
mresln <- MISO(gene, reads=reads, noIterations=3000, noBurnIn=2500,
               noLag=10, start="linear", noChains=5)

cor(exp, postMean(mres))
cor(exp, postMean(mres2))
cor(exp, postMean(mresln))
sum(abs(exp-postMean(mres)))
sum(abs(exp-postMean(mres2)))
sum(abs(exp-postMean(mresln)))
percentAccepted(mres)
percentAccepted(mres2)
percentAccepted(mresln)

#################

gene3 <- selectGenes(mus, 3)

exp3 <- norm((3:1)^2)
reads3 <- simulateReads(gene3, readLength=35, noReads=4000, expression=exp3)
rexp3 <- realPsi(gene3, reads=reads3)

mres3 <- MISO(gene3, reads=reads3, noIterations=5000,
              noBurnIn=2500, noLag=10)
mres32 <- MISO(gene3, reads=reads3, noChains=5, noIterations=2000,
              noBurnIn=1500, noLag=10)
mresln3 <- MISO(gene3, reads=reads3, noIterations=5000,
              noBurnIn=2500, noLag=10, start="linear")

cor(exp3, postMean(mres3))
cor(exp3, postMean(mres32))
cor(exp3, postMean(mresln3))
sum(abs(exp3-postMean(mres3)))
sum(abs(exp3-postMean(mres32)))
sum(abs(exp3-postMean(mresln3)))
percentAccepted(mres3)
percentAccepted(mres32)
percentAccepted(mresln3)

#########

gene2 <- selectGenes(mus, 23)

exp2 <- norm((2:1)^2)
reads2 <- simulateReads(gene2, readLength=35, noReads=4000, expression=exp2)
rexp2 <- realPsi(gene2, reads=reads2)

mresx2 <- MISO(gene2, reads=reads2, noIterations=5000,
               noBurnIn=2500, noLag=10)
mresx22 <- MISO(gene2, reads=reads2, noChains=5, noIterations=2000,
                noBurnIn=1500, noLag=100)
mresln2 <- MISO(gene2, reads=reads2, noIterations=5000,
                noBurnIn=2500, noLag=10, start="linear")

cor(exp2, postMean(mresx2))
cor(exp2, postMean(mresx22))
cor(exp2, postMean(mresln2))
sum(abs(exp2-postMean(mresx2)))
sum(abs(exp2-postMean(mresx22)))
sum(abs(exp2-postMean(mresln2)))
percentAccepted(mresx2)
percentAccepted(mresx22)
percentAccepted(mresln2)
