
library(splicing)

gffFile <- system.file("test_files/Atp2b1.mm9.gff", package="splicing")

gff <- readGFF3(gffFile)

.Call("R_splicing_noexons_one", gff, 1L, PACKAGE="splicing")

