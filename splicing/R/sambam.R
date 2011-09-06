
## TODO: attributes

readSAM <- function(filename, region=NULL) {
  if (is.null(region)) {
    .Call("R_splicing_read_sambam", as.character(filename),
          PACKAGE="splicing")
  } else {
    .Call("R_splicing_read_sambam_region", as.character(filename),
          as.character(region), PACKAGE="splicing")
  }
}

.SQ <- list('Mus musculus'='@SQ\tSN:mm9_allJxns.1e6.withSG.35.4.fa\tLN:547718976
@SQ\tSN:chrRibo\tLN:45309
@SQ\tSN:chr1\tLN:197195432
@SQ\tSN:chr2\tLN:181748087
@SQ\tSN:chr3\tLN:159599783
@SQ\tSN:chr4\tLN:155630120
@SQ\tSN:chr5\tLN:152537259
@SQ\tSN:chr6\tLN:149517037
@SQ\tSN:chr7\tLN:152524553
@SQ\tSN:chr8\tLN:131738871
@SQ\tSN:chr9\tLN:124076172
@SQ\tSN:chrM\tLN:16299
@SQ\tSN:chrX\tLN:166650296
@SQ\tSN:chr10\tLN:129993255
@SQ\tSN:chr11\tLN:121843856
@SQ\tSN:chr12\tLN:121257530
@SQ\tSN:chr13\tLN:120284312
@SQ\tSN:chr14\tLN:125194864
@SQ\tSN:chr15\tLN:103494974
@SQ\tSN:chr16\tLN:98319150
@SQ\tSN:chr17\tLN:95272651
@SQ\tSN:chr18\tLN:90772031
@SQ\tSN:chr19\tLN:61342430')

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


  attr(tab, "header") <- splicing.i.SAMheader(reads)
  tab
}

writeSAM <- function(reads, conn)
  UseMethod("writeSAM")

## TODO: attributes

writeSAM.splicingSAM <- function(reads, conn) {
  header <- splicing.i.SAMheader(reads)
  tab <- sprintf("%s\t%i\t%s\t%i\t%i\t%s\t%s\t%s\t%i\t%s\t%s",
                 as.character(reads$qname), as.integer(reads$flag),
                 as.character(reads$chrname[reads$chr+1]),
                 as.integer(reads$position), as.integer(reads$mapq),
                 as.character(reads$cigar), as.character(reads$rnext),
                 as.character(reads$pairpos), as.integer(reads$tlen),
                 as.character(reads$seq), as.character(reads$qual))
  cat(file=conn, header, "\n", sep="", paste(tab, collapse="\n"), "\n")
}
