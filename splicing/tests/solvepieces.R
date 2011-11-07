
library(splicing)

gene <- createGene(list(c(1,20), c(31,40)),
                   list(c(1), c(1,2)))

exons <- .Call("R_splicing_gff_exon_start_end", gene, 1L, PACKAGE="splicing")

.Call("R_splicing_numeric_cigar", exons, noIso(gene), 1L)

exons$start <- exons$start - 1L
exons$end   <- exons$end   - 1L

.Call("R_splicing_numeric_cigar", exons, noIso(gene), 0L)
