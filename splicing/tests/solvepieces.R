
library(splicing)

gene <- createGene(list(c(1,20), c(31,40)),
                   list(c(1), c(1,2)))

exons <- .Call("R_splicing_gff_exon_start_end", gene, 1L, PACKAGE="splicing")

.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 0L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 1L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 19L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 20L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 21L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 25L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 29L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 30L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 31L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 40L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L, 50L)

exons$start <- exons$start - 1L
exons$end   <- exons$end   - 1L

.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 0L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 1L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 19L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 20L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 21L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 25L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 29L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 30L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 31L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 40L)
.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L, 50L)
