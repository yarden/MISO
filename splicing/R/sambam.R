
isReads <- function(reads) {
  inherits(reads, "splicingSAM")
}

readSAM <- function(filename, region=NULL, geneStructure=NULL, gene=1) {
  if (!is.null(region) && !is.null(geneStructure)) {
    stop("At most one of `region' and `geneStructure' can be given")
  }
  if (is.null(region) && is.null(geneStructure)) {
    .Call("R_splicing_read_sambam", as.character(filename),
          PACKAGE="splicing")
  } else if (!is.null(region)) {
    .Call("R_splicing_read_sambam_region", as.character(filename),
          as.character(region), PACKAGE="splicing")
  } else {
    region <- getRegion(geneStructure, gene)
    .Call("R_splicing_read_sambam_region", as.character(filename),
          as.character(region), PACKAGE="splicing")
  }
}

SAMFile2BAMFile <- function(samfile, bamfile) {
  res <- .Call("R_splicing_sam2bam", as.character(samfile),
               as.character(bamfile),
               PACKAGE="splicing")
  invisible(res)
}

sortBAMFile <- function(file, outprefix, key=c("position", "qname")) {
  key <- match.arg(key)
  key <- switch(key, "position"=0L, "qname"=1L)
  res <- .Call("R_splicing_bam_sort", as.character(file),
               as.character(outprefix), key, PACKAGE="splicing")
  invisible(res)
}

indexBAMFile <- function(file) {
  res <- .Call("R_splicing_bam_index", as.character(file),
               PACKAGE="splicing")
  invisible(res)
}
                        
noReads <- function(reads)
  UseMethod("noReads")

noReads.splicingSAM <- function(reads) {
  length(reads$position)
}

print.splicingSAM <- function(x, ...) {
  reads <- x
  if (isPaired(reads)) {
    mess <- sprintf("SAM/BAM %i read pairs, %i single reads, %i sequences",
                    noPairs(reads), noSingles(reads),
                    length(seqNames(reads)))
  } else {
    mess <- sprintf("SAM/BAM %i single-end reads, %i sequences",
                    noSingles(reads), length(seqNames(reads)))
  }
  cat(mess, "\n")
}

isPaired <- function(reads)
  UseMethod("isPaired")

isPaired.splicingSAM <- function(reads) {
  reads$paired
}

isMixed <- function(reads)
  UseMethod("isMixed")

isMixed.splicingSAM <- function(reads) {
  reads$noPairs > 0 && reads$noSingles > 0
}

noPairs <- function(reads)
  UseMethod("noPairs")

noPairs.splicingSAM <- function(reads) {
  reads$noPairs
}

noSingles <- function(reads)
  UseMethod("noSingles")

noSingles.splicingSAM <- function(reads) {
  reads$noSingles
}

seqNames <- function(reads)
  UseMethod("seqNames")

seqNames.splicingSAM <- function(reads) {
  reads$chrname
}

eachSeqName <- function(reads)
  UseMethod("eachSeqName")

eachSeqName.splicingSAM <- function(reads) {
  reads$chrname[reads$chr+1]
}

getReadLength <- function(reads)
  UseMethod("getReadLength")

getReadLength.splicingSAM <- function(reads) {
  len <- strsplit(reads$cigar, "[A-Z]")
  typ <- lapply(strsplit(reads$cigar, "[0-9]+"), "[", -1)
  res <- mapply(len, typ, FUN=function(l, t) sum(as.numeric(l[t=="M"])))
  sort(unique(res))
}

## To support different read lengths within the same file

## TODO: speed this up

eachReadLength <- function(reads) {
  len <- strsplit(reads$cigar, "[A-Z]")
  typ <- lapply(strsplit(reads$cigar, "[0-9]+"), "[", -1)
  res <- mapply(len, typ, FUN=function(l, t) sum(as.numeric(l[t=="M"])))
}  

toTable <- function(reads, ...)
  UseMethod("toTable")

splicing.i.SAMheader <- function(reads) {
  paste(sep="\n", collapse="\n",
        "@HD\tVN:1.0\tSO:unsorted",
        sprintf("@SQ\tSN:%s\tLN:%i", reads$chrname, reads$chrlen))
} 

## TODO: attributes

toTable.splicingSAM <- function(reads, ...) {
  attrnames <- listAttributes(reads)
  if (length(attrnames) == 0) { 
    tab <- data.frame(QNAME=reads$qname,
                      FLAG=reads$flag,
                      RNAME=reads$chrname[reads$chr+1],
                      POS=reads$position,
                      MAPQ=reads$mapq,
                      CIGAR=reads$cigar,
                      RNEXT=reads$rnext,
                      PNEXT=reads$pairpos,
                      TLEN=reads$tlen,
                      SEQ=reads$seq,
                      QUAL=reads$qual,
                      ...)
  } else {
    tab <- data.frame(QNAME=reads$qname,
                      FLAG=reads$flag,
                      RNAME=reads$chrname[reads$chr+1],
                      POS=reads$position,
                      MAPQ=reads$mapq,
                      CIGAR=reads$cigar,
                      RNEXT=reads$rnext,
                      PNEXT=reads$pairpos,
                      TLEN=reads$tlen,
                      SEQ=reads$seq,
                      QUAL=reads$qual,
                      ATTR=getAttributes(reads),
                      ...)
  }    


  attr(tab, "header") <- splicing.i.SAMheader(reads)
  tab
}

writeSAM <- function(reads, conn)
  UseMethod("writeSAM")

writeSAM.splicingSAM <- function(reads, conn) {
  header <- splicing.i.SAMheader(reads)
  attr <- getAttributes(reads)
  if (is.null(attr)) { attr <- "" } else { attr <- paste("\t", sep="", attr) }
  tab <- sprintf("%s\t%i\t%s\t%i\t%i\t%s\t%s\t%s\t%i\t%s\t%s%s",
                 as.character(reads$qname), as.integer(reads$flag),
                 as.character(reads$chrname[reads$chr+1]),
                 as.integer(reads$position), as.integer(reads$mapq),
                 as.character(reads$cigar), as.character(reads$rnext),
                 as.character(reads$pairpos), as.integer(reads$tlen),
                 as.character(reads$seq), as.character(reads$qual),
                 attr)
  cat(file=conn, header, "\n", sep="", paste(tab, collapse="\n"), "\n")
}

delbit <- function(v, bit) {
  .Call("R_splicing_delbit", v, bit,
        PACKAGE="splicing")
}

getbit <- function(v, bit) {
  .Call("R_splicing_getbit", v, bit,
        PACKAGE="splicing")
}

getStrand <- function(reads)
  UseMethod("getStrand")

getStrand.splicingSAM <- function(reads) {
  getbit(reads$flag, 5L)
}

getAttributes <- function(reads)
  UseMethod("getAttributes")

getAttributes.splicingSAM <- function(reads) {
  reads$attributes
}

listAttributes <- function(reads)
  UseMethod("listAttributes")

listAttributes.splicingSAM <- function(reads) {
  if (is.null(reads$attributes)) {
    character()
  } else { 
    a <- strsplit(reads$attributes, "\t", fixed=TRUE)
    unique(unlist(lapply(a, function(x) substr(x, 1, 2))))
  }
}

getAttribute <- function(reads, attr)
  UseMethod("getAttribute")

getAttribute.splicingSAM <- function(reads, attr) {
  a <- reads$attributes
  if (is.null(a)) { return(a) }
  a <- strsplit(a, "\t", fixed=TRUE)
  attr <- paste("^", attr, ":", sep="")
  sapply(a, function(x) {
    g <- grep(attr, x)
    if (length(g)==0) {
      NA
    } else {
      sp <- strsplit(x[g[1]], ":", fixed=TRUE)[[1]]
      if (sp[2] %in% c("i", "f")) {
        as.numeric(sp[3])
      } else {
        sp[3]
      }
    }
  })
}

getIsoform <- function(reads)
  UseMethod("getIsoform")

getIsoform.splicingSAM <- function(reads) {
  getAttribute(reads, attr="XI")
}

selectReads <- function(reads, idx)
  UseMethod("selectReads")

selectReads.splicingSAM <- function(reads, idx) {
  no <- noReads(reads)
  idx <- sort(seq_len(no)[idx])
  if (!is.null(attr <- reads$attributes)) { attr <- attr[idx] }
  res <- list(chrname=reads$chrname, chrlen=reads$chrlen,
              chr=reads$chr[idx], qname=reads$qname[idx],
              cigar=reads$cigar[idx], position=reads$position[idx],
              flag=reads$flag[idx], pairpos=reads$pairpos[idx],
              noPairs=0L, noSingles=0L, paired=reads$paired,
              mapq=reads$mapq[idx], rnext=reads$rnext[idx],
              tlen=reads$tlen[idx], seq=reads$seq[idx], qual=reads$qual[idx],
              mypair=reads$mypair[idx], attributes=attr)
  no2 <- length(res$cigar)
  class(res) <- class(reads)

  ## Fix 1) number of pairs, 2) number of singles, 3) paired, 4) flags,
  ## 5) pairpos and 6) mate pointers

  if (isPaired(reads)) {
    tran <- rep(-1L, no)
    tran[idx] <- seq_along(idx)-1L

    res$mypair <- sapply(reads$mypair[idx],
                         function(x) if (x==-1) { -1 } else { tran[x+1] })
    res$noPairs <- as.integer(sum(res$mypair != -1)/2)
    res$noSingles <- as.integer(no2 - 2 * res$noPairs)
    res$paired <- res$noPairs != 0
    res$pairpos[res$mypair==-1] <- 0
    res$flag[res$mypair==-1] <- delbit(res$flag[res$mypair==-1], 1L)
    
  } else {
    res$noSingles <- as.integer(no2)
  }
  
  res
}

filterReads <- function(reads, ...)
  UseMethod("filterReads")

## Remove all paired-end reads

filter.noPaired <- function(reads) {
  selectReads(reads, which(reads$pairpos==0))
}

## Remove all single-end reads

filter.noSingle <- function(reads) {
  selectReads(reads, which(reads$pairpos!=0))
  reads
}

## Remove paired-end reads where the mates are
## not on the opposite strand. Based on attributes.

filter.notOppositeStrand <- function(reads) {
  ## Keep single reads
  keep1 <- which(reads$pairpos==0)

  ## And the appropriate paired reads
  paired <- which(reads$pairpos!=0)
  strand <- getStrand(reads)
  keep2 <- paired[which(strand[paired] != strand[reads$mypair[paired]+1])]
  
  selectReads(reads, sort(c(keep1, keep2)))
}

## Remove reads that do not match any isoform of a given
## gene

filter.notMatching <- function(reads, gff, paired=isPaired(reads), ...) {
  
  mat <- matchIso(geneStructure=gff, reads=reads, paired=paired, ...)

  if (paired) { 
    keep <- which(colSums(mat[[1]]) != 0)
  } else {
    keep <- which(colSums(mat) != 0)
  }
  
  selectReads(reads, keep)
}

## Remove inconsistent reads, where the flags do not match
## the rest of the read

filter.badFlags <- function(reads) {

  pflag <- getbit(reads$flag, 1L)
  keep <- which((!pflag & reads$pairpos==0) |
                ( pflag & reads$pairpos!=0))
  
  selectReads(reads, keep)
}

## Keep only a given region, on a given chromosome

filter.keepRegion <- function(reads, sequence, start, end) {
  rl <- eachReadLength(reads)
  keep <- which(reads$chrname[reads$chr+1] == sequence &
                reads$position >= start-rl+1 & reads$position <= end)
  selectReads(reads, keep)
}

## Filter out too short reads

filter.dropShort <- function(reads, lengthLimit) {
  rl <- eachReadLength(reads)
  keep <- which(rl >= lengthLimit)
  selectReads(reads, keep)
}

filter.dropLong <- function(reads, lengthLimit) {
  rl <- eachReadLength(reads)
  keep <- which(rl <= lengthLimit)
  selectReads(reads, keep)
}  

filter.dropOtherLength <- function(reads, length) {
  rl <- eachReadLength(reads)
  keep <- which(rl == length)
  selectReads(reads, keep)
}

myfilters <- list(noPaired=filter.noPaired,
                  noSingle=filter.noSingle,
                  notOppositeStrand=filter.notOppositeStrand,
                  notMatching=filter.notMatching,
                  badFlags=filter.badFlags,
                  keepRegion=filter.keepRegion,
                  dropShort=filter.dropShort,
                  dropLong=filter.dropLong,
                  dropOtherLength=filter.dropOtherLength)

FILTER <- function(name, ...) {
  if (! name %in% names(myfilters)) {
    stop("Unknown filter: ", name)
  }
  res <- list(name=name, args=list(...))
  class(res) <- "splicingReadFilter"
  res
}

filterReads.splicingSAM <- function(reads, ..., verbose=TRUE) {

  is.filter <- function(x) {
    inherits(x, "splicingReadFilter")
  }
  
  filters <- list(...)
  if ( !all(sapply(filters, is.character) | sapply(filters, is.filter)) ) {
    stop("Only character and FILTER() objects are allowed as filters")
  }

  sapply(filters[sapply(filters, is.character)], function(x) {
    if (length(x)!=1) { stop("Invalid filter: ", x) }
    if (!x %in% names(myfilters)) { stop("Unknown filter: ", x) }
  })
  
  for (i in seq_along(filters)) {
    origlen <- noReads(reads)
    if (is.character(filters[[i]])) {
      name <- filters[[i]]
      reads <- myfilters[[ filters[[i]] ]](reads)
    } else {
      name <- filters[[i]]$name
      reads <- do.call(myfilters[[ filters[[i]]$name ]],
                       c(list(reads=reads), filters[[i]]$args))
    }
    newlen <- noReads(reads)
    if (verbose) {
      message(sprintf("`%s' filtered out %i reads, %i remaining",
                      name, origlen-newlen, newlen))
    }
  }

  reads
}

mergeReads <- function(...) {
  reads <- list(...)
  chrname <- unique(unlist(lapply(reads, "[[", "chrname")))

  allchrname <- unlist(lapply(reads, "[[", "chrname"))
  allchrlen <- unlist(lapply(reads, "[[", "chrlen"))
  chrlen <- allchrlen[match(chrname, allchrname)]

  chr <- match(unlist(lapply(reads, function(x) x$chrname[x$chr+1])),
               chrname)-1
  qname <- unlist(lapply(reads, "[[", "qname"))
  if (any(duplicated(unlist(lapply(reads, function(x) unique(x$qname)))))) {
    warning("Duplicate QNAME fields while merging reads")
  }
  cigar <- unlist(lapply(reads, "[[", "cigar"))
  position <- unlist(lapply(reads, "[[", "position"))
  flag <- unlist(lapply(reads, "[[", "flag"))
  pairpos <- unlist(lapply(reads, "[[", "pairpos"))
  noPairs <- sum(sapply(reads, "[[", "noPairs"))
  noSingles <- sum(sapply(reads, "[[", "noSingles"))
  paired <- any(sapply(reads, "[[", "paired"))
  mapq <- unlist(lapply(reads, "[[", "mapq"))
  rnext <- unlist(lapply(reads, "[[", "rnext"))
  tlen <- unlist(lapply(reads, "[[", "tlen"))
  seq <- unlist(lapply(reads, "[[", "seq"))
  qual <- unlist(lapply(reads, "[[", "qual"))

  offset <- c(0, cumsum(sapply(reads, noReads))[-length(reads)])
  mypair <- unlist(mapply(reads, offset, SIMPLIFY=FALSE,
                          FUN=function(r, o) {
                            ifelse(r$mypair==-1, -1, r$mypair + o)
                          })
                   )

  if (any(sapply(reads, function(x) !is.null(x$attributes)))) {
    attributes <- unlist(lapply(reads, function(x) {
      if (is.null(x$attributes)) character(noReads(x)) else x$attributes
    }))
  } else {
    attributes <- NULL
  }
  sampleProb <- NULL

  res <- list(chrname=chrname, chrlen=chrlen, chr=chr, qname=qname,
              cigar=cigar, position=position, flag=flag, pairpos=pairpos,
              noPairs=noPairs, noSingles=noSingles, paired=paired,
              mapq=mapq, rnext=rnext, tlen=tlen, seq=seq, qual=qual,
              mypair=mypair, attributes=attributes, sampleProb=sampleProb)
  class(res) <- "splicingSAM"
  res
}
                  
estimateFragLength <- function(genes, readsfile, reads, min_length, ...) {
  if (isGFF3(genes)) {
    genes <- constitutiveExons(genes, min_length, ...)
  }
  if (!inherits(genes, "splicingExonset")) {
    stop("`genes' must be an exon set or a GFF3 object")
  }
  if (missing(readsfile)) { readsfile <- NULL }
  if (missing(reads))     { reads     <- NULL }

  if (is.null(readsfile) && is.null(reads)) {
    stop("Please give the reads or the name of the reads file.")
  }
  if (!is.null(readsfile) && !is.null(reads)) {
    stop("Please give only one of `readsfile' and `reads'")
  }
  
  .Call("R_splicing_estimate_fraglength", genes, readsfile, reads,
        PACKAGE="splicing")
}
