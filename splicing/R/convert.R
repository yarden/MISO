
listEnsemblGenes <- function(baseURL="ftp://ftp.ensembl.org/pub/current_gtf/") {

  require(RCurl)
  
  list <- getURL(baseURL)
  list <- strsplit(list, "\n")[[1]]
  list <- read.table(textConnection(list), header=FALSE)
  species <- as.character(list[[ncol(list)]])
  species <- paste(toupper(substring(species,1,1)), substring(species,2),
                   sep="")
  species <- gsub("_", " ", species)
  species
}

## Function to download a GTF file from ensemble and
## convert it to a GFF3.

getEnsemblGenes <- function(speciesName, 
                            baseURL="ftp://ftp.ensembl.org/pub/current_gtf/",
                            verbose=TRUE) {

  require(RCurl)

  ## Create the URL
  
  speciesName2 <- gsub(" ", "_", tolower(speciesName))
  list <- getURL(baseURL)
  if (!any(grepl(speciesName2, list))) {
    stop("Unknown species")
  }
  
  list2 <- getURL(paste(baseURL, sep="/", speciesName2, ""))
  fname <- tail(strsplit(list2, " ")[[1]], 1)

  ## Download the file
  
  tmp <- tempfile()
  if (verbose) { message("Downloading file") }
  download.file(paste(baseURL, sep="/", speciesName2, fname), tmp,
                quiet=!verbose)
  if (verbose) { message("Uncompresssing and parsing file") }
  ## TODO: something faster
  myfile <- gzfile(tmp)
  lines <- read.delim(myfile, header=FALSE, stringsAsFactors=FALSE)
  unlink(tmp)
  
  cn <- c("seqname", "source", "feature", "start", "end", "score", "strand",
          "frame")
  if (ncol(lines) < 8) {
    warning("Less than eight columns in GTF file")
    cn <- cn[seq_len(ncol(lines))]
  }
  if (ncol(lines) > 8) { cn <- c(cn, "attributes") }
  if (ncol(lines) > 9) { cn <- c(cn, "comments") }
  if (ncol(lines) > 10) { warning("More than ten columns in GTF file") }
  colnames(lines) <- cn
  attr(lines, "species") <- speciesName
  
  gtf2gff3(lines)
}

gtf2gff3 <- function(gtf, verbose=TRUE) {

  ## Empty input is special case
  if (nrow(gtf)==0) {
    res <- list(seqid=character(), source=character(), type=character(),
                start=integer(), end=integer(), score=character(),
                strand=character(), phase=character(),
                attributes=character())
    class(res) <- c("gff3", "data.frame")
    return(res)
  }

  gtf$seqname    <- as.character(gtf$seqname)
  gtf$source     <- as.character(gtf$source)
  gtf$feature    <- as.character(gtf$feature)
  gtf$score[gtf$score=="."] <- "-1"
  gtf$score      <- as.numeric(gtf$score)
  gtf$strand     <- as.character(gtf$strand)  
  gtf$frame      <- as.character(gtf$frame)
  gtf$attributes <- as.character(gtf$attributes)
  
  gid <- sub(".*gene_id ([^;]+);?.*", "\\1", gtf$attributes)
  tid <- sub(".*transcript_id ([^;]+);?.*", "\\1", gtf$attributes)
  ord <- order(gtf$seqname, gid, tid, gtf$feature)
  gtf <- gtf[ord,]

  res <- .Call("R_splicing_gtf2gff3", gtf, gid[ord], tid[ord],
               PACKAGE="splicing")
  res <- res[ res$start != -1, ]
  
  class(res) <- c("gff3")

  types <- c(gene=SPLICING_GENE, mRNA=SPLICING_MRNA,
             exon=SPLICING_EXON, CDS=SPLICING_CDS,
             start_codon=SPLICING_START_CODON,
             stop_codon=SPLICING_STOP_CODON)
  strands <- c("+"=SPLICING_STRAND_PLUS,
               "-"=SPLICING_STRAND_MINUS,
               "."=SPLICING_STRAND_UNKNOWN)
             
  res$type <- types[res$type]
  res$phase[res$phase=="."] <- "-1"
  res$phase <- as.integer(res$phase)
  res$species <- attr(gtf, "species")
  
  res <- addGidTid(res)

  res$seqid_str <- unique(res$seqid)
  res$seqid <- match(res$seqid[res$gid+1L], res$seqid_str)-1L
  res$source_str <- unique(res$source)
  res$source <- match(res$source[res$gid+1L], res$source_str)-1L
  res$strand <- strands[res$strand[res$gid+1L]]

  res
}

addGidTid <- function(gff3) {
  gff3$gid <- which(gff3$type==SPLICING_GENE)-1L
  gff3$tid <- which(gff3$type==SPLICING_MRNA)-1L
  if (! "ID" %in% names(gff3)) {
    gff3$ID <- sub(".*ID=([^;]+);?.*", "\\1", gff3$attributes)
  }
  if (! "parent" %in% names(gff3)) {
    gff3$parent <- sub(".*;Parent=([^;]+);?.*", "\\1", gff3$attributes)
    gff3$parent <- match(gff3$parent, gff3$ID)-1L
    gff3$parent[is.na(gff3$parent)] <- -1L
  }
  gff3
}

readGFF3 <- function(file) {
  .Call("R_splicing_read_gff", as.character(file),
        PACKAGE="splicing")
}

writeGFF3 <- function(gff3, file) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  res <- .Call("R_splicing_write_gff", gff3, file,
               PACKAGE="splicing")
  invisible(res)
}

## Gives number of genes

noGenes <- function(gff3)
  UseMethod("noGenes")

noGenes.gff3 <- function(gff3) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  length(gff3$gid)
}

## Sequence ids

seqIds <- function(gff3)
  UseMethod("seqIds")

seqIds.gff3 <- function(gff3) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  res <- gff3$seqid_str[gff3$seqid+1]
  names(res) <- geneIds(gff3)
  res
}  

## Gives the gene ids

geneIds <- function(gff3)
  UseMethod("geneIds")

geneIds.gff3 <- function(gff3) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  gff3$ID[gff3$gid+1]
}

## Select some genes, based on index, or ID

selectGenes <- function(gff3, idx)
  UseMethod("selectGenes")

selectGenes.gff3 <- function(gff3, idx) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }

  if (is.character(idx)) {
    idx <- match(idx, gff3$ID)
    if (any(is.na(idx))) {
      stop("unknown gene selected")
    }
  }
  
  start <- gff3$gid+1
  end <- c(start[-1]-1, length(gff3$start))
  sel <- unlist(mapply(start[idx], end[idx], FUN=":", SIMPLIFY=FALSE))

  res <- with(gff3, {
    list(seqid_str=seqid_str, seqid=seqid[idx],
         source_str=source_str, source=source[idx],
         type=type[sel], start=start[sel], end=end[sel],
         score=score[sel], strand=strand[idx], phase=phase[sel],
         attributes=attributes[sel], ID=ID[sel], parent=parent[sel])
  })

  gp <- res$parent >= 0
  res$parent[gp] <- as.integer(match(gff3$ID[res$parent[gp]+1], res$ID)-1)

  class(res) <- "gff3"
  res <- addGidTid(res)
  res
}

## Number of isoforms

noIso <- function(gff3)
  UseMethod("noIso")

noIso.gff3 <- function(gff3) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  .Call("R_splicing_gff_noiso", gff3, PACKAGE="splicing")
}

selectIso <- function(gff3, iso)
  UseMethod("selectIso")

selectIso.gff3 <- function(gff3, iso) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  if (noGenes(gff3) != 1) {
    stop("GFF should contain a single gene")
  }

  ktid <- gff3$tid[iso]
  ktidID <- gff3$ID[ktid+1L]
  par <- sapply(gff3$parent, function(x) if (x==-1) NA else gff3$ID[x+1])
  keep <- which(gff3$type==SPLICING_GENE | gff3$ID %in% ktidID |
                par %in% ktidID)

  attributes <- if (is.null(gff3$attributes)) NULL else gff3$attributes[keep]
  
  res <- list(seqid_str=gff3$seqid_str, seqid=gff3$seqid,
              source_str=gff3$source_str, source=gff3$source,
              type=gff3$type[keep], start=gff3$start[keep],
              end=gff3$end[keep], score=gff3$score[keep],
              strand=gff3$strand, phase=gff3$phase[keep],
              attributes=attributes, ID=gff3$ID[keep])
  res$parent <- match(par[keep], res$ID)-1L
  res$parent[is.na(res$parent)] <- -1L
  res$gid <- which(res$type == SPLICING_GENE)-1L
  res$tid <- which(res$type == SPLICING_MRNA)-1L
  class(res) <- "gff3"
  res
}  

## Number of exons, exons are different if they have a different
## start position

noExons <- function(gff3)
  UseMethod("noExons")

noExons.gff3 <- function(gff3) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  gs <- gff3$gid+1
  ex <- which(gff3$type == SPLICING_EXON)
  sign <- paste(sep="-", gff3$start[ex], gff3$end[ex])
  ex <- ex[ !duplicated(sign) ]
  ggr <- cut(ex, breaks=c(gs, length(gff3$start)+1))
  as.numeric(unname(table(ggr)))
}

## Total length of all exons. Overlapping portions are counted
## multiple times.

totalExonLength <- function(gff3, overlap=TRUE)
  UseMethod("totalExonLength")

## TODO: rewrite in C to make it faster, especially the overlap=FALSE
##       case

totalExonLength.gff3 <- function(gff3, overlap=TRUE) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  gs <- gff3$gid
  ex <- which(gff3$type == SPLICING_EXON)
  if (overlap) {
    si <- paste(sep="-", gff3$start[ex], gff3$end[ex])
    un <- which(!duplicated(si))
    le <- gff3$end[ex][un] - gff3$start[ex][un] + 1
    gr <- cut(ex[un], breaks=c(gs, length(gff3$start)+1))
    res <- tapply(le, gr, sum)
  } else {
    gr <- cut(ex, breaks=c(gs, length(gff3$start)+1))
    res <- tapply(ex, gr, function(eidx) {
      myex <- data.frame(start=gff3$start[eidx],
                         end=gff3$end[eidx])
      myex <- myex[order(myex$start, myex$end),,drop=FALSE]
      len <- 0
      aex <- c(myex$start[1], myex$end[1])
      for (i in seq_len(nrow(myex))[-1]) {
        nex <- c(myex$start[i], myex$end[i])
        if (nex[1] > aex[2]) {
          len <- len + aex[2]-aex[1]+1
          aex <- nex
        } else {
          aex[2] <- nex[2]
        }
      }
      len <- len + aex[2]-aex[1]+1
      len
    })
  }
  names(res) <- geneIds(gff3)
  res
}

## Get the names of isoforms for all genes

getIso <- function(gff3, collapse=FALSE)
  UseMethod("getIso")

getIso.gff3 <- function(gff3, collapse=FALSE) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  g <- gff3$gid+1
  mr <- gff3$tid+1
  tid <- gff3$ID[mr]
  ggr <- cut(mr, breaks=c(g, length(gff3$start)+1))
  res <- tapply(tid, ggr, c, simplify=FALSE)
  if (collapse) {
    res <- unname(unlist(res))
  } else {
    names(res) <- gff3$ID[g]
  }
  res
}

## Is this a GFF3 object?

isGFF3 <- function(object) {
  inherits(object, "gff3")
}

## What type of genes?

geneTypes <- function(gff3)
  UseMethod("geneTypes")

geneTypes.gff3 <- function(gff3) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  res <- gff3$source_str[gff3$source+1]
  names(res) <- geneIds(gff3)
  res
}

## The length of the gene(s)

geneLength <- function(gff3)
  UseMethod("geneLength")

geneLength.gff3 <- function(gff3) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  g <- gff3$gid+1
  res <- gff3$end[g] - gff3$start[g] + 1
  names(res) <- geneIds(gff3)
  res
}

## The length of the different isoforms, after splicing

isoLength <- function(gff3)
  UseMethod("isoLength")

isoLength.gff3 <- function(gff3) {
  .Call("R_splicing_gff_isolength", gff3, PACKAGE="splicing")
}

getSpecies <- function(gff3)
  UseMethod("getSpecies")

getSpecies.gff3 <- function(gff3) {
  res <- gff3$species
  if (is.null(res)) NA else res
}

constitutiveExons <- function(gff3, min.length, mode=c("full", "all"))
  UseMethod("constitutiveExons")

constitutiveExons.gff3 <- function(gff3, min.length, mode=c("full", "all")) {
  if (!isGFF3(gff3)) {
    stop("Not a GFF3 object")
  }
  mode <- match.arg(mode)
  mode <- switch(mode, "full"=0L, "all"=1L)
  .Call("R_splicing_constitutive_exons", gff3, as.integer(min.length),
        as.integer(mode), PACKAGE="splicing")
}

getExonStart <- function(gff3, gene=1) {
  mygff <- selectGenes(gff3, gene[1])
  mr <- mygff$tid+1
  ex <- which(mygff$type==SPLICING_EXON)
  res <- tapply(ex, cut(ex, breaks=c(mr, length(mygff$start)+1)),
                function(x) mygff$start[x], simplify=FALSE)
  names(res) <- getIso(mygff)[[1]]
  res
}

## TODO: make it faster

getExonEnd <- function(gff3, gene=1) {
  mygff <- selectGenes(gff3, gene[1])
  mr <- mygff$tid
  ex <- which(mygff$type==SPLICING_EXON)
  res <- tapply(ex, cut(ex, breaks=c(mr, length(mygff$start)+1)),
                function(x) mygff$end[x], simplify=FALSE)
  names(res) <- getIso(mygff)[[1]]
  res
}

## TODO: make it faster

getExonLength <- function(gff3, gene=1) {
  mygff <- selectGenes(gff3, gene[1])
  mr <- gff3$tid+1
  ex <- which(mygff$type==SPLICING_EXON)
  res <- tapply(ex, cut(ex, breaks=c(mr, length(mygff$start)+1)),
                function(x) mygff$end[x]-mygff$start[x]+1, simplify=FALSE)
  names(res) <- getIso(mygff)[[1]]
  res
}

getExonSig <- function(gff3, gene=1)
  UseMethod("getExonSig")

getExonSig.gff3 <- function(gff3, gene=1) {
  st <- getExonStart(gff3, gene)
  en <- getExonEnd(gff3, gene)
  mapply(st, en, FUN=function(a,b) paste(a, b, sep=":", collapse=";"))
}

getRegion <- function(gff3, gene=1)
  UseMethod("getRegion")

getRegion.gff3 <- function(gff3, gene=1) {
  if (length(gene) != 1) {
    stop("`gene' must be a scalar")
  }
  start <- gff3$start[gff3$gid[gene]+1]
  end <- gff3$end[gff3$gid[gene]+1]
  seqid <- gff3$seqid[gene]
  seqid_str <- gff3$seqid_str[seqid+1]
  paste(sep="", seqid_str, ":", start, "-", end)  
}

print.gff3 <- function(x, verbose=TRUE, ...) {
  gff3 <- x
  nog <- noGenes(gff3)
  spec <- getSpecies(gff3)
  if (is.na(spec)) { spec <- "Unknown species" }
  if (nog != 1 || !verbose) {
    cat(sprintf('GFF3 %s, %i genes, %i transcripts.\n', spec, nog,
                length(gff3$tid)))
  } else {
    cat(sprintf('GFF3 %s gene, %i isoforms.\n', spec, length(gff3$tid)))
    ## TODO: more
  }
}

print.splicingExonset <- function(x, ...) {
  s <- if (length(x$start)==1) "" else "s"
  cat(sprintf('Set of %i exon%s, from %i sequences.\n', length(x$start), s,
              length(x$seqid_str)))
  invisible(x)
}

##############
'
GFF3 Caenorhabditis_elegans, 45461 genes, 54959 transcripts.
'
