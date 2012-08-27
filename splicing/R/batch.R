
run <- function(runFunction, writeFunction, geneStructure, readsfile, 
                results = c("return", "files", "Rfiles"),
                resultDir = ".", overWrite = FALSE, readLength = NULL,
                verbose = TRUE, snowCluster = NULL,
                seed = NULL, dropBadCigar=FALSE, ...) {

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
    myclusterExport(snowCluster, c("readsfile",
                                   "readLength", "verbose", "ids",
                                   "geneStructure", "results", "overWrite"))
    myapply <- function(...) parLapply(cl=snowCluster, ...)
  } else {
    myapply <- lapply
  }

  res <- myapply(seq_along(ids), function(x) {

    if (results=="files") {
      fname <- paste(sep="", resultDir, "/", ids[x], ".miso")
    } else if (results=="Rfiles") {
      fname <- paste(sep="", resultDir, "/", ids[x], ".Rdata")            
    }
    
                 
    if (!overWrite && results %in% c("files", "Rfiles") &&
        file.exists(fname)) {
      return(fname)
    }
    
    region <- getRegion(geneStructure, x)
    reads <- readSAM(readsfile, region=region)
    if (dropBadCigar) {
      reads <- selectReads(reads, grep("[^MNSHDI\\=X0-9]", reads$cigar,
                                       invert=TRUE))
    }
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
      writeFunction(res, fname)
      fname
    } else if (results=="Rfiles") {
      save(res, file=fname)
      fname
    }
  })
  
  names(res) <- ids
  res
}
