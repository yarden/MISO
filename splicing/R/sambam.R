
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

writeSAM <- function(sam, file, SQ, organism, ...) {

  if (!inherits(sam, "sam")) {
    stop("Not a SAM object")
  }

  if (is.character(file)) {
    file <- file(file, "wt")
    toClose <- TRUE
  } else {
    toClose <- FALSE
  }

  ## Header
  HD <- sprintf("@HD\tVN:%s\tSO:%s", attr(sam, "HD_VN"),
                attr(sam, "HD_SO"))
  PG <- sprintf("@PG\tID:%s\tVN:%s\tCL:%s", attr(sam, "PG_ID"),
                attr(sam, "PG_VN"), attr(sam, "PG_CL"))
  cat(HD, file=file, sep="\n")

  ## SQ lines, this is organism specific information
  if (!missing(SQ)) {
    cat(SQ, file=file, sep="\n")
  } else if (!missing(organism)) {
    if (! organism %in% names(.SQ)) {
      stop("Unknown organism")
    }
    cat(.SQ[[organism]], file=file, sep="\n")
  } else {
    warning("@SQ headers are not written, give `SQ' or `organism'")
  }

  cat(PG, file=file, sep="\n")
  
  ## Data
  if (attr(sam, "paired")) {
    QNAME <- rep(paste(attr(sam, "QNAME"), sep=".", seq_len(nrow(sam)/2)),
                 each=2)
  } else {
    QNAME <- paste(attr(sam, "QNAME"), sep=".", seq_len(nrow(sam)))
  }
  tab <- data.frame(QNAME=QNAME,
                    FLAG=sam$FLAG,
                    RNAME=sam$RNAME,
                    POS=sam$position,
                    MAPQ=255,
                    CIGAR=sam$cigar,
                    RNEXT='=',
                    PNEXT=if ('PNEXT' %in% names(sam)) sam$PNEXT else 0,
                    TLEN= if ('TLEN' %in% names(sam)) sam$TLEN else 0,
                    SEQ='*',
                    QUAL='*',
                    XI=paste(sep="", "XI:i:", sam$isoform, "\tXS:A:",
                      ifelse(bitAnd(16, sam$FLAG), "-", "+")))

  write.table(tab, file=file, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE, ...)

  if (toClose) { close(file) }

  invisible(NULL)
}

sam2bam <- function(samfile, bamdir, misoDir,
                    python="/usr/bin/env python") {

  cmd <- sprintf("%s '%s/sam_to_bam.py' --convert %s %s",
                 python, misoDir, samfile, bamdir)
  system(cmd)  
}
