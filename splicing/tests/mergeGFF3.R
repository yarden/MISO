
library(splicing)

gffFile1 <- system.file("test_files/Atp2b1.mm9.gff", package="splicing")
gffFile2 <- system.file("test_files/ENSMUSG00000026173.gff",
                        package="splicing")

g1 <- readGFF3(gffFile1)
g2 <- readGFF3(gffFile2)

g <- merge(g1, g2)
str(g)

