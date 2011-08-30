
plotIso <- function(geneStructure, gene=geneIds(geneStructure)[1],
                    xlab="", ylab="", ...) {

  geneStructure <- selectGenes(geneStructure, gene)

  ## The empty plot
  plot(NA, type="n", xlim=c(0,1), ylim=c(0,1), axes=FALSE, frame=FALSE,
       xlab=xlab, ylab=ylab, ...)

  ## How many isoforms, how long the gene
  lines <- noIso(geneStructure)[1]
  len <- geneLength(geneStructure)[1]
  
  ## Labels
  tid <- getIso(geneStructure)[[1]]
  labwidth <- max(strwidth(paste(tid, " ")))
  text(x=0, y=seq(1,by=2,length=lines)/2/lines,
       adj=c(0,1/2), rev(tid))

  mr <- which(geneStructure$type==SPLICING_MRNA)
  ex <- which(geneStructure$type==SPLICING_EXON)
  st <- which(geneStructure$type==SPLICING_START_CODON)
  en <- which(geneStructure$type==SPLICING_STOP_CODON)

  wex <- tapply(ex, cut(ex, breaks=c(mr, length(geneStructure$start)+1)), c)
  wst <- tapply(st, cut(st, breaks=c(mr, length(geneStructure$start)+1)), c)
  wen <- tapply(en, cut(en, breaks=c(mr, length(geneStructure$start)+1)), c)

  .scale <- function(val, fx1, fx2, tx1, tx2) {
    (val-fx1) / (fx2-fx1) * (tx2-tx1) + tx1
  }

  .droplast <- function(x) x[-length(x)]
  
  .plot <- function(lowerleft, upperright, exons, startcodon, stopcodon,
                    startpos, endpos, strand) {

    xsize <- upperright[1]-lowerleft[1]
    ysize <- upperright[2]-lowerleft[2]

    exons <- exons[ order(geneStructure$start[exons]) ]
    
    estart <- geneStructure$start[exons]
    eend <- geneStructure$end[exons]
    
    start <- end <- 0
    if (length(startcodon)>0 && length(stopcodon)>0 &&
        !is.na(startcodon[1]) && !is.na(stopcodon[1])) {
      if (strand=="-") {
        start <- geneStructure$start[stopcodon[1]]
        end <- geneStructure$end[startcodon[1]]
      } else {
        start <- geneStructure$start[startcodon[1]]
        end <- geneStructure$end[stopcodon[1]]
      }
    }
    
    ## Plot exons first
    x <- .scale((eend+estart)/2, startpos, endpos,
                lowerleft[1], upperright[1])
    rw <- (eend-estart) / (endpos-startpos) * xsize
    rh <- ysize/4
    symbols(inches=FALSE, add=TRUE,
            x=x, y=rep((upperright[2]+lowerleft[2])/2, length(x)),
            rectangles=cbind(rw,rh), fg="black", bg="black")

    ## Transcribed exons
    for (e in seq_along(exons)) {
      if (start <= estart[e] && eend[e] <= end) {
        ## whole exon
        rw <- (eend[e]-estart[e]) / (endpos-startpos) * xsize
        rh <- ysize/2
        symbols(inches=FALSE, add=TRUE,
                x=.scale((eend[e]+estart[e])/2, startpos, endpos,
                  lowerleft[1], upperright[1]),
                y=(upperright[2]+lowerleft[2])/2,
                rectangles=cbind(rw, rh), fg="black", bg="lightgrey")
      } else if (estart[e] <= start && start <= eend[e] && eend[e] <= end) {
        ## end of exon
        rw <- (eend[e]-start) / (endpos-startpos) * xsize
        rh <- ysize/2
        symbols(inches=FALSE, add=TRUE,
                x=.scale((eend[e]+start)/2, startpos, endpos,
                  lowerleft[1], upperright[1]),
                y=(upperright[2]+lowerleft[2])/2,
                rectangles=cbind(rw, rh), fg="black", bg="lightgrey")
      } else if (start <= estart[e] && estart[e] <= end && end <= eend[e]) {
        ## start of exon
        rw <- (end-estart[e]) / (endpos-startpos) * xsize
        rh <- ysize/2
        symbols(inches=FALSE, add=TRUE,
                x=.scale((end+estart[e])/2, startpos, endpos,
                  lowerleft[1], upperright[1]),
                y=(upperright[2]+lowerleft[2])/2,
                rectangles=cbind(rw, rh), fg="black", bg="lightgrey")        
      } else if (estart[e] <= start && end <= eend[e]) {
        ## middle of exon
        rw <- (end-start) / (endpos-startpos) * xsize
        rh <- ysize/2
        symbols(inches=FALSE, add=TRUE,
                x=.scale((end+start)/2, startpos, endpos,
                  lowerleft[1], upperright[1]),
                y=(upperright[2]+lowerleft[2])/2,
                rectangles=cbind(rw, rh), fg="black", bg="lightgrey")        
      }
    }
    
    ## Connecting lines
    if (length(exons) > 1) {
      x0 <- .scale(eend[-length(exons)],
                   startpos, endpos, lowerleft[1], upperright[1])
      y0 <- ifelse(eend > start & eend <= end, 6/8, 5/8) * ysize +
        lowerleft[2]
      y02 <- ifelse(estart >= start & estart < end, 6/8, 5/8) * ysize +
        lowerleft[2]
      x1 <- .scale(geneStructure$start[exons][-1],
                   startpos, endpos, lowerleft[1], upperright[1])
      y1 <- rep(lowerleft[2]+ysize*7/8)
      segments(x0=x0, y0=.droplast(y0), x1=(x0+x1)/2, y1=y1)
      segments(x0=(x0+x1)/2, y0=y1, x1=x1, y1=y02[-1])
    }
  }
  
  for (i in 1:lines) {
    .plot(lowerleft=c(labwidth, (i-1)/lines),
          upperright=c(1, i/lines),
          exons=wex[[lines-i+1]], startcodon=wst[[lines-i+1]],
          stopcodon=wen[[lines-i+1]],
          startpos=geneStructure$start[1], endpos=geneStructure$end[1],
          strand=geneStructure$strand[1])
  }

  invisible(geneStructure)
}

plotIsoSize <- function(geneStructure, gene=geneIds(geneStructure)[1],
                        labels=NULL, stripwidth=5, stripheight=1/5) {

  geneStructure <- selectGenes(geneStructure, gene)

  no <- noIso(geneStructure)[1]
  if (is.null(labels)) {
    labels <- getIso(geneStructure)[[1]]
  }
  labwidth <- max(strwidth(paste(labels, " "), units="inches"))
  width <- labwidth + stripwidth
  height <- stripheight * no
  c(width, height)
}

plotIsoPDF <- function(geneStructure, file, gene=geneIds(geneStructure)[1],
                       labelwidth=NULL, mar=c(0,0,2,0), ...) {

  pdf(tmp <- tempfile())
  plot.new()
  size <- plotIsoSize(geneStructure, gene, labels=labelwidth)
  dev.off()
  unlink(tmp)
  
  pdf(file, width=size[1], height=size[2])
  par(mar=mar)
  plotIso(geneStructure, gene, ...)
  dev.off()
}
