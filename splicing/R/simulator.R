
## Return the start positions of the exons, for all genes
## TODO: make it faster

SPLICING_GENE <- 0L
SPLICING_MRNA <- 1L
SPLICING_EXON <- 2L
SPLICING_START_CODON <- 3L
SPLICING_STOP_CODON <- 4L

SPLICING_STRAND_PLUS <- 0L
SPLICING_STRAND_MINUS <- 1L
SPLICING_STRAND_UNKNOWN <- 2L

getExonStart <- function(gff3, gene) {
  mygff <- selectGenes(gff3, gene[1])
  mr <- mygff$tid+1
  ex <- which(mygff$type==SPLICING_EXON)
  res <- tapply(ex, cut(ex, breaks=c(mr, length(mygff$start)+1)),
                function(x) mygff$start[x], simplify=FALSE)
  names(res) <- getIso(mygff)[[1]]
  res
}

## TODO: make it faster

getExonEnd <- function(gff3, gene) {
  mygff <- selectGenes(gff3, gene[1])
  mr <- mygff$tid
  ex <- which(mygff$type==SPLICING_EXON)
  res <- tapply(ex, cut(ex, breaks=c(mr, length(mygff$start)+1)),
                function(x) mygff$end[x], simplify=FALSE)
  names(res) <- getIso(mygff)[[1]]
  res
}

## TODO: make it faster

getExonLength <- function(gff3, gene) {
  mygff <- selectGenes(gff3, gene[1])
  mr <- gff3$tid+1
  ex <- which(mygff$type==SPLICING_EXON)
  res <- tapply(ex, cut(ex, breaks=c(mr, length(mygff$start)+1)),
                function(x) mygff$end[x]-mygff$start[x]+1, simplify=FALSE)
  names(res) <- getIso(mygff)[[1]]
  res
}

## Generate reads for a single gene

simulateReads <- function(geneStructure, gene=1, expression,
                          noReads, readLength, paired=FALSE, fragmentProb=NULL,
                          fragmentStart=0L, normalMean, normalVar, numDevs) {

  if (!paired) {
    res <- .Call("R_splicing_simulate_reads", geneStructure,
                 as.integer(gene), as.double(expression),
                 as.integer(noReads), as.integer(readLength),
                 PACKAGE="splicing")
    res$paired <- FALSE
  } else {
    if (!is.null(fragmentProb)) { fragmentProb <- as.double(fragmentProb) }
    res <- .Call("R_splicing_simulate_paired_reads", geneStructure,
                 as.integer(gene), as.double(expression),
                 as.integer(noReads), as.integer(readLength),
                 fragmentProb, as.integer(fragmentStart), as.double(normalMean),
                 as.double(normalVar), as.double(numDevs))
    res$paired <- TRUE
  }

  res
}

## TODO: implement 'geneStructure'

createGene <- function(exons, isoforms, id="insilicogene",
                       seqid="seq1", source="protein_coding",
                       strand=c("+", "-", "."), geneStructure=NULL) {

  exons <- lapply(exons, as.integer)
  isoforms <- lapply(isoforms, as.integer)
  strand <- match.arg(strand)
  strand <- switch(strand, "+"=0L, "-"=1L, "."=2L)
  .Call("R_splicing_create_gene", exons, isoforms, id, seqid, source,
        strand,
        PACKAGE="splicing")
}

mergeGenes <- function(..., mode=c("disjunct")) {
  genes <- list(...)
  if (length(genes)==0) {
    stop("No genes are given")
  }
  if (any(!sapply(genes, isGFF3))) {
    stop("Genes must be GFF3 objects")
  }
  allgids <- unlist(sapply(genes, geneIds))
  if (any(duplicated(allgids))) {
    stop("Duplicate gene ids: ", allgids[duplicated(allgids)][1])
  }
  alltids <- unlist(sapply(genes, getIso))
  if (any(duplicated(alltids))) {
    stop("Duplicate transcript ids: ", alltids[duplicated(alltids)][1])
  }
  mode <- match.arg(mode)
  if (length(genes)==1) {
    return(genes[[1]])
  }

  if (mode=="disjunct") {
    maxpos <- sapply(genes, function(x) max(x$end))
    maxpos <- cumsum(maxpos)-maxpos[1]
    for (i in seq_along(genes)) {
      genes[[i]]$start <- genes[[i]]$start + maxpos[i]
      genes[[i]]$end   <- genes[[i]]$end   + maxpos[i]
    }
    res <- do.call(rbind, genes)
  }

  res <- addGidTid(res)
  res
}

## Reads per kilobase. If 'overlap' is true, then the overlapping
## parts of overlapping exons are counted multiple times. Otherwise
## the length of the mRNA is calculated by considering each base that
## is included in some mature isoform.

rpk <- function(genes, reads, overlap=TRUE) {

  if (length(unique(seqIds(genes))) != 1 &&
      ! 'RNAME' %in% names(reads)) {
    stop("Genes from multiple sequences, but no sequence IDs in reads")
  }
  
  exlen <- totalExonLength(genes, overlap=overlap)
  
  myrpk <- function(g) {
    mygene <- selectGenes(genes, g)
    myreads <- reads[mygene$start[1] <= reads$position  &
                     reads$position  <= mygene$end[1],,drop=FALSE]
    mat <- matchIso(mygene, reads=myreads)
    myreads <- myreads[!grepl("^0+$", mat),,drop=FALSE]
    nrow(myreads) / exlen[g] * 1000
  }

  sapply(seq_len(noGenes(genes)), myrpk)
}

## Calculate the real expression profile. This might differ a bit from
## the one that was prescribed for the read generator

realPsi <- function(gene, reads, readLength=33) {
  isolen <- isoLength(gene)[[1]]
  effLen <- isolen - readLength + 1
  ii <- as.character(1:noIso(gene)[1])
  it <- table(reads$isoform)[ii]
  it[is.na(it)] <- 0
  res <- as.vector(it / effLen)
  res / sum(res)
}
