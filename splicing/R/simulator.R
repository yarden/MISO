
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
                          noReads, readLength, paired=FALSE, insertProb=NULL,
                          insertStart=0L, normalMean, normalVar, numDevs) {

  if (!paired) {
    res <- .Call("R_splicing_simulate_reads", geneStructure,
                 as.integer(gene), as.double(expression),
                 as.integer(noReads), as.integer(readLength),
                 PACKAGE="splicing")
    res$paired <- FALSE
  } else {
    if (!is.null(insertProb)) { insertProb <- as.double(insertProb) }
    res <- .Call("R_splicing_simulate_paired_reads", geneStructure,
                 as.integer(gene), as.double(expression),
                 as.integer(noReads), as.integer(readLength),
                 insertProb, as.integer(insertStart), as.double(normalMean),
                 as.double(normalVar), as.double(numDevs))
    res$paired <- TRUE
  }

  res
}

genReadsForGene <- function(geneStructure, expression,
                            gene=geneIds(geneStructure)[1],
                            noReads=NA, coverage=10, readLength=33,
                            paired=FALSE, insertDistr=NULL,
                            normalInsert=c(250,15)) {
  
  geneStructure <- selectGenes(geneStructure, gene)
  no <- noIso(geneStructure)  
  
  if (noGenes(geneStructure) != 1) {
    stop("`geneStructure' must have exactly one gene")
  }
  if (no == 0) {
    stop("No isoforms in `geneStructure'")
  }
  if (noExons(geneStructure) == 0) {
    stop("No exons in `geneStructure'")
  }

  if (length(expression) != no) {
    stop("Length of `expression' must match number of isoforms")
  }

  if (paired && is.null(insertDistr) && is.null(normalInsert)) {
    stop("Insert length distribution is not given for paired end reads")
  }
  
  ## Calculate the insert length distribution, if not explicitly given
  ## TODO: use 2 sd rather
  if (paired  && is.null(insertDistr) && !is.null(normalInsert)) {
    insertDistr <- diff(pnorm(0:(2*normalInsert[1]), mean=normalInsert[1],
                              sd=normalInsert[2]))
  }
  
  ## Calculate the number of samples to generate. For this we use
  ## the expression levels of the isoforms, and their length as well.

  isolen <- isoLength(geneStructure)[[1]]
  expression <- expression / sum(expression)

  if (is.na(noReads)) {
    meanLen <- sum(expression * isolen)
    if (!paired) {
      noSamples <- ceiling(meanLen * coverage / readLength)
    } else {
      noSamples <- ceiling(meanLen * coverage / (readLength*2))
    }
  } else {
    noSamples <- noReads
  }
  
  ## Sampling, which isoform the read is coming from. Probability is
  ## proportional to the length of the isoform and to its expression
  ## level
  if (!paired) {
    sampleProb <- ifelse(isolen >= readLength,
                         (isolen-readLength+1) * expression, 0)
  } else {
    cinp <- cumsum(insertDistr)
    startp <- sapply(isolen, function(L) {
      if (length(cinp) < L-2*readLength) {
        tp <- c(rep(1, L-2*readLength-length(cinp)), rev(cinp))
      } else {
        tp <- rev(cinp[1:(L-2*readLength)])
      }
      sum(tp)
    })
    sampleProb <- ifelse(isolen >= readLength*2, startp * expression, 0)
  }
  whichIso <- sample(seq_along(isolen), replace=TRUE, noSamples,
                     prob=sampleProb)

  ## Place the genes on the original sequence, sort the exons first
  exstart <- getExonStart(geneStructure, gene=1)
  exend   <- getExonEnd(geneStructure, gene=1)
  exlen   <- getExonLength(geneStructure, gene=1)
  exord <- lapply(exstart, order)
  exstart <- mapply(exstart, exord, FUN=function(x, o) x[o], SIMPLIFY=FALSE)
  exend   <- mapply(exend, exord, FUN=function(x, o) x[o], SIMPLIFY=FALSE)
  exlen   <- mapply(exlen, exord, FUN=function(x, o) x[o], SIMPLIFY=FALSE)

  .droplast <- function(x) x[-length(x)]
  
  shift <- mapply(exstart, exend, SIMPLIFY=FALSE, FUN=function(es, ee) {
    cumsum(as.numeric(es)) - c(0,.droplast(cumsum(as.numeric(ee)))) -
      0:(length(es)-1)-1
  })

  exlim <- lapply(exlen, function(x) c(1, cumsum(x)+1))
  
  ## Decide where the read starts within the isoform, for single end
  ## reads this is easy, because it is uniform.
  ## For paired end reads we need to take into account the
  ## insert length distribution. In this case we store the actual
  ## insert length as well.
  if (!paired) {
    readStart <- sapply(whichIso, function(w) {
      spos <- sample(1:(isolen[w]-readLength+1), 1)
      ex <- which(exlim[[w]] > spos)[1]-1
      spos + shift[[w]][ex]
    })
  } else {
    readStart <- sapply(whichIso, function(w) {
      posRL <- 1:(isolen[w]-readLength*2)
      if (length(posRL) > length(cinp)) {
        prob <- c(rep(1, length(posRL)-length(cinp)), rev(cinp))
      } else {
        prob <- rev(cinp[1:length(posRL)])
      }
      spos <- sample(posRL, 1, prob=prob)
      ex <- which(exlim[[w]] > spos)[1]-1
      pil <- 1:(isolen[w]-spos+1-readLength*2)
      prob2 <- insertDistr[pil]
      prob2 <- ifelse(is.na(prob2), 0, prob2)
      il <- sample(pil, 1, prob=prob2)
      spos2 <- spos+readLength+il
      ex2 <- which(exlim[[w]] > spos2)[1]-1
      ## if (is.na(spos2 + shift[[w]][ex2])) { browser() }
      c(spos + shift[[w]][ex], spos2 + shift[[w]][ex2], il)
    })
    ## Trick, first n elements will be the normal read start coordinates
    iLen <- readStart[3,]
    readStart <- t(readStart[1:2,])
    
  }

  ## Need to calculate the CIGAR string, to show the skipped and
  ## matched parts of the reads
  cigar <- sapply(seq_along(whichIso), function(i) {
    w <- whichIso[i]
    ex <- which(exend[[w]] >= readStart[i])[1]
    rs <- readStart[i]
    rl <- readLength
    cig <- ""
    while (exend[[w]][ex] < rs+rl-1) {
      cig <- paste(cig, sep="", exend[[w]][ex]-rs+1, "M",
                   exstart[[w]][ex+1]-exend[[w]][ex]-1, "N")
      rl <- rl - (exend[[w]][ex] - rs + 1)
      rs <- exstart[[w]][ex+1]
      ex <- ex + 1
    }
    cig <- paste(cig, sep="", rl, "M")
    cig
  })
  if (paired) {
    cigar2 <- sapply(seq_along(whichIso), function(i) {
      w <- whichIso[i]
      ex <- which(exend[[w]] >= readStart[i,2])[1]
      rs <- readStart[i,2]
      rl <- readLength
      cig <- ""
      while (exend[[w]][ex] < rs+rl-1) {
        cig <- paste(cig, sep="", exend[[w]][ex]-rs+1, "M",
                     exstart[[w]][ex+1]-exend[[w]][ex]-1, "N")
        rl <- rl - (exend[[w]][ex] - rs + 1)
        rs <- exstart[[w]][ex+1]
        ex <- ex + 1
      }
      cig <- paste(cig, sep="", rl, "M")
      cig
    })
  }

  if (!paired) {
    res <- data.frame(isoform=whichIso, position=as.integer(readStart),
                      cigar=cigar,
                      FLAG=if (geneStructure$strand[1]=="+") 0 else 16,
                      RNAME=geneStructure$seqid[1],
                      stringsAsFactors=FALSE)
  } else {
    FLAG <- if (geneStructure$strand[1]=="+") c(67,131) else c(83,147)
    position <- as.vector(t(readStart))
    PNEXT <- as.vector(t(readStart[,2:1]))
    res <- data.frame(isoform=rep(whichIso, each=2),
                      position=as.integer(position),
                      cigar=as.vector(rbind(cigar,cigar2)),
                      FLAG=rep(FLAG, length(whichIso)),
                      RNAME=geneStructure$seqid[1],
                      PNEXT=PNEXT,
                      TLEN=ifelse(PNEXT > position,
                        PNEXT-position+readLength,
                        PNEXT-position-readLength),
                      stringsAsFactors=FALSE)
  }
  class(res) <- c("sam", "data.frame")

  attr(res, "geneId") <- geneIds(geneStructure)
  attr(res, "HD_VN") <- "1.0"
  attr(res, "HD_SO") <- "unsorted"
  attr(res, "PG_ID") <- "splicing R package simulator"
  attr(res, "PG_VN") <- packageVersion("splicing")
  attr(res, "PG_CL") <- "TODO"
  attr(res, "QNAME") <- "TODO"
  attr(res, "FLAG")  <- if (geneStructure$strand[1]=="+") 0 else 16
  attr(res, "paired") <- paired
    
  res
}

genReads <- function(geneStructure, expression, genes=geneIds(geneStructure),
                     ...) {

  if (!isGFF3(geneStructure)) {
    stop("Not a GFF3 object")
  }

  if (!is.list(expression)) {
    stop("`expression' should be a ist")
  }

  geneStructure <- selectGenes(geneStructure, genes)
  
  if (noGenes(geneStructure) != length(expression)) {
    stop("`expression' length should match length of `genes'")
  }

  if (any(duplicated(genes))) {
    warning("Duplicated gene ids in `genes'")
  }

  if (!is.null(names(expression))) {
    if (any(duplicated(names(expression)))) {
      warning("Duplicated gene ids in `expression'")
    }
  } else {
    names(expression) <- geneIds(geneStructure)
  }

  if (any(sort(names(expression)) != sort(geneIds(geneStructure)))) {
    stop("Gene ids in `expression' and `genes'")
  }
      
  res <- lapply(geneIds(geneStructure), function(id) {
    tmp <- genReadsForGene(geneStructure, expression[[id]], gene=id, ...)
    tmp$geneId <- attr(tmp, "geneId")
    if (!"FLAG" %in% names(tmp)) { tmp$FLAG <- attr(tmp, "FLAG") }
    if (!"PNEXT" %in% names(tmp)) { tmp$PNEXT <- 0 }
    attr(tmp, "geneId") <- attr(tmp, "FLAG") <- NULL
    tmp
  })
    
  res <- do.call(rbind, res)

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
