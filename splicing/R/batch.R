
run <- function(runFunction, writeFunction, geneStructure, readsfile, 
                results = c("return", "files", "Rfiles"),
                resultDir = ".", overWrite = FALSE, readLength = NULL,
                verbose = TRUE, snowCluster = NULL,
                seed = NULL, ...) {

  results <- match.arg(results)

  if (results=="files" && is.null(writeFunction)) {
    stop("`results=\"files\"' is not implemented for this solver")
  }
  
  if (is.character(geneStructure)) {
    geneStructure <- readGFF3(geneStructure)
  }
  if (!isGFF3(geneStructure)) {
    stop("Invalid gene structure(s), must be a GFF3 object or a filename")
  }

  ids <- geneIds(geneStructure)
  
  get.region <- function(gene) {
    start <- geneStructure$start[geneStructure$gid[gene]+1]
    end <- geneStructure$end[geneStructure$gid[gene]+1]
    seqid <- geneStructure$seqid[gene]
    seqid_str <- geneStructure$seqid_str[seqid+1]
    paste(sep="", seqid_str, ":", start, "-", end)
  }
  
  myclusterExport <- function(cl, list) {
    for (name in list) {
      clusterCall(cl, assign, name, get(name))
    }
  }

  if (!is.null(snowCluster)) {
    require(snow)
    if (!is.null(seed)) {
      clusterCall(snowCluster, set.seed, seed)
    }
    myclusterExport(snowCluster, c("get.region", "readsfile",
                                   "readLength", "verbose", "ids",
                                   "geneStructure", "results", "overWrite"))
    myapply <- function(...) parLapply(cl=snowCluster, ...)
  } else {
    myapply <- lapply
  }

  res <- myapply(seq_along(ids), function(x) {

    if (results=="files") {
      fname <- paste(sep="", resultDir, "/", ids[gene], ".miso")
    } else if (results=="Rfiles") {
      fname <- paste(sep="", resultDir, "/", ids[gene], ".Rdata")            
    }
    
                 
    if (!overWrite && results %in% c("files", "Rfiles") &&
        file.exists(fname)) {
      return(fname)
    }
    
    region <- get.region(x)
    reads <- readSAM(readsfile, region=region)
    rl <- if (is.null(readLength)) {
      getReadLength(reads)
    } else {
      readLength
    }
    if (verbose) { message("Running gene # ", x, ", ", ids[x]) }
    res <- runFunction(geneStructure, gene=x, reads=reads, readLength=rl, ...)

    if (results=="return") {
      res
    } else if (results=="files") {
      writeFunction(misoResult, fname)
      fname
    } else if (results=="Rfiles") {
      save(misoResult, file=fname)
      fname
    }
  })
  
  names(res) <- ids
  res
}
