
library(splicing)

gffFile <- system.file("test_files/Atp2b1.mm9.gff", package="splicing")

gff <- readGFF3(gffFile)

.Call("R_splicing_noexons_one", gff, 1L, PACKAGE="splicing")

########

library(splicing)

gene <- createGene(list(c(11,30), c(41,50)),
                   list(c(1), c(1,2)))

.Call("R_splicing_genomic_to_iso_all", gff=gene, gene=1L, position=11L,
      PACKAGE="splicing")
# 1 1
.Call("R_splicing_genomic_to_iso_all", gff=gene, gene=1L, position=41L,
      PACKAGE="splicing")
# -1 21
.Call("R_splicing_genomic_to_iso_all", gff=gene, gene=1L, position=35L,
      PACKAGE="splicing")
# -1 -1
.Call("R_splicing_genomic_to_iso_all", gff=gene, gene=1L, position=51L,
      PACKAGE="splicing")
# -1 -1
.Call("R_splicing_genomic_to_iso_all", gff=gene, gene=1L, position=5L,
      PACKAGE="splicing")
# -1 -1

.Call("R_splicing_iso_to_genomic_all", gff=gene, gene=1L, position=1L,
      PACKAGE="splicing")
# 11 11
.Call("R_splicing_iso_to_genomic_all", gff=gene, gene=1L, position=20L,
      PACKAGE="splicing")
# 30 30
.Call("R_splicing_iso_to_genomic_all", gff=gene, gene=1L, position=30L,
      PACKAGE="splicing")
# -1 50
.Call("R_splicing_iso_to_genomic_all", gff=gene, gene=1L, position=31L,
      PACKAGE="splicing")
# -1 -1
