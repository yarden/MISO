
readSAM <- function(filename, region=NULL) {
  if (is.null(region)) {
    .Call("R_splicing_read_sambam", as.character(filename),
          PACKAGE="splicing")
  } else {
    .Call("R_splicing_read_sambam_region", as.character(filename),
          as.character(region), PACKAGE="splicing")
  }
}

noReads <- function(reads)
  UseMethod("noReads")

noReads.splicingSAM <- function(reads) {
  length(reads$position)
}

print.splicingSAM <- function(reads) {
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

readLength <- function(reads)
  UseMethod("readLength")

readLength.splicingSAM <- function(reads) {
  len <- strsplit(reads$cigar, "[A-Z]")
  typ <- lapply(strsplit(reads$cigar, "[0-9]+"), "[", -1)
  res <- mapply(len, typ, FUN=function(l, t) sum(as.numeric(l[t=="M"])))
  sort(unique(res))
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
  res <- list(chrname=reads$chrname, chrlen=reads$chrlen,
              chr=reads$chr[idx], qname=reads$chr[idx],
              cigar=reads$cigar[idx], position=reads$position[idx],
              flag=reads$flag[idx], pairpos=reads$pairpos[idx],
              noPairs=0L, noSingles=0L, paired=reads$paired,
              mapq=reads$mapq[idx], rnext=reads$rnext[idx],
              tlen=reads$tlen[idx], seq=reads$seq[idx], qual=reads$qual[idx],
              mypair=reads$mypair[idx])
  no2 <- length(res$cigar)
  class(res) <- class(reads)

  ## Fix 1) number of pairs, 2) number of singles, 3) paired, 4) flags,
  ## 5) pairpos and 6) mate pointers

  if (isPaired(reads)) {
    keep <- sort(seq_len(no)[idx])
    tran <- integer(no)
    tran[keep] <- seq_along(keep)
    tran <- tran-1L

    res$mypair <- tran[reads$mypair+1]
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

filter.notMatching <- function(reads, gff, paired=reads$paired, ...) {
  
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

myfilters <- list(noPaired=filter.noPaired,
                  noSingle=filter.noSingle,
                  notOppositeStrand=filter.notOppositeStrand,
                  notMatching=filter.notMatching,
                  badFlags=filter.badFlags)

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
