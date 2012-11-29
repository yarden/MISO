
plotIso <- function(geneStructure, gene=1, xlab="", ylab="", ...) {

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
      if (strand==SPLICING_STRAND_MINUS) {
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

plotIsoSize <- function(geneStructure, gene=1,
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

plotIsoPDF <- function(geneStructure, file, gene=1,
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

## TODO:
##  - stacked area charts

plotReadsSize <- function(gene, reads, misoResult=NULL) {
  width   <- 7 + if (is.null(misoResult)) 0 else 2
  height  <- .15 + length(reads) * .75 + .35 + noIso(gene) * .5 + .5
  c(width=width, height=height)
}

plotReads <- function(gene, reads, misoResult=NULL,
                      ## Colors
                      sampleColors="#e78800ff",
                      isoformColors=rainbow_hcl(noIso(gene)),
                      ## Default sizes of various plot regions
                      titleHeight=0.15, axisHeight=0.35,
                      geneBaseHeight=0.5, geneIsoHeight=0.5,
                      misoWidth=2) {

  require(igraph)
  require(colorspace)
  
  ## Check arguments
  if (!isGFF3(gene)) { stop("`gene' must be a GFF3 object") }
  if (any(!sapply(reads, isReads))) {
    stop("reads must be a list of reads")
  }
  if (!is.null(misoResult)) {
    if (length(misoResult) != length(reads)) {
      stop("Number of samples and length of MISO results must match")
    }
  }

  noSamples <- length(reads)

  sampleColors <- rep(sampleColors, length.out=noSamples)
  isoformColors <- rep(isoformColors, length.out=noIso(gene))

  ## Layout
  laymat <- cbind(seq_len(noSamples+3))
  laymat <- cbind(laymat, laymat+max(laymat)-1)
  laymat[1,2] <- 1

  din <- par("din")
  layh <- c(title=titleHeight, rep(NA, noSamples), axis=axisHeight,
            geneStructure=(noIso(gene)*geneIsoHeight + geneBaseHeight))
  histh <- (din[2] - sum(layh, na.rm=TRUE)) / noSamples
  layh[ is.na(layh) ] <- histh

  mw <- if (is.null(misoResult)) 0 else misoWidth
  layw <- c(din[1]-mw, mw)

  if (layw[1] <= 0) {
    stop("Figure width too small, minimum suggested width is 3 inches")
  }
  if (layw[1] <= 1.5) {
    warning("Figure width too small, minimum suggested width is 3 inches")
  }

  if (any(layh <= 0)) {
    stop("Figure height too small, minimum suggested height for this\n",
         "gene and samples is ", .5 * noSamples+.15+.35+noIso(gene)*.5+.5,
         " inches")
  }
  if (histh < .5) {
    warning("Figure height too small, minimum suggested height for this ",
            "gene and samples is ", .5 * noSamples+.15+.35+noIso(gene)*.5+.5,
            " inches")
  }
  layout(laymat, widths=layw, heights=layh)

  ## Title
  par(mar=c(0,0,0,0))
  plot.new()
  text(sum(par("usr")[1:2])/2, sum(par("usr")[3:4])/2, cex=1, geneIds(gene),
       adj=c(1/2,1), xpd=NA)
  print("title:"); print(par("fin"))
  
  ## Common parameters
  start <- getExonStart(gene)
  end <- getExonEnd(gene)
  xlim <- range(unlist(start), unlist(end))
  
  ## Create histogram for reads
  ## TODO: write a member function for splicingSAM instead of this
  ## We also count the number of reads across junctions
  rhist <- lapply(reads, function(r) {
    no <- rep(0, xlim[2]-xlim[1]+1)
    junc <- matrix(nrow=0, ncol=2)
    for (i in seq_along(r$pos)) {
      pos <- r$pos[i]
      cigar <- r$cigar[i]
      len <- as.numeric(strsplit(cigar, "[A-Z]")[[1]])
      typ <- lapply(strsplit(cigar, "[0-9]+"), "[", -1)[[1]]
      idx <- pos - xlim[1]+1
      for (j in seq_along(len)) {
        if (typ[j]=="M") {
          idx2 <- idx:(idx+len[j]-1)
          no[idx2] <- no[idx2] + 1
          idx <- idx + len[j]
        } else if (typ[j]=="N") {
          junc <- rbind(junc, c(idx-1, idx+len[j]))
          idx <- idx + len[j]
        } else {
          stop("Unknown CIGAR string character")
        }
      }
    }
    list(no=no, junc=junc)
  })

  junc <- lapply(rhist, "[[", "junc")
  rhist <-lapply(rhist, "[[", "no")

  ## Omit the zero values, plot it as several separate polygons
  mypolygon <- function(x, y, ...) {
    idx <- 1
    for (i in seq_along(x)) {
      if ((y[i] == 0 || i==length(x)) && i > 1  && y[i-1] != 0) {
        polygon(c(x[idx], x[idx:i], x[i]), c(0, y[idx:i], 0), ...)
        idx <- i
      } else if (y[i]==0 && i<length(x) && y[i+1] != 0) {
        idx <- i+1
      }
    }
  }

  ## Where to put the junctions
  jpos <- lapply(junc, function(jj) {
    uj <- unique(jj)
    el <- c()
    if (nrow(uj)==0) { return(numeric()) }
    if (nrow(uj)==1) { return(1) }
    for (i in 1:(nrow(uj)-1)) {
      for (j in (i+1):nrow(uj)) {
        if ( ! (uj[i,1] > uj[j,2] || uj[j,1] > uj[i,2]) ) {
          el <- c(el, i, j)
        }
      }
    }
    g <- graph(el, n=nrow(uj), directed=FALSE)
    bm <- bipartite.mapping(g)
    if (!bm$res) {
      warning("Cannot plot all junctions")
      jp <- rep(0:1, length.out=nrow(uj))
    } else {
      jp <- as.numeric(bm$type)
    }
    if (sum(jp==0) > sum(jp==1)) { jp <- 1-jp }
    jp
  })
  
  ## Plot the samples
  ylim <- c(0, max(unlist(rhist))*1.3)  
  for (i in seq_along(reads)) {
    par(mar=c(0,5,1,1)+.1)
    plot(NA, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="# reads",
         axes=FALSE)
    xval <- seq(min(unlist(start)), max(unlist(end)))
    mypolygon(xval, rhist[[i]], col=sampleColors[i], border=NA)
    axis(2, las=1, cex.axis=.8, lwd=.4)
    allj <- unique(junc[[i]])
    jc <- table(match(paste(junc[[i]][,1], sep=":", junc[[i]][,2]),
                      paste(allj[,1], sep=":", allj[,2])))
    for (j in seq_len(nrow(allj))) {
      pos <- jpos[[i]][j]
      jfromx <- allj[j,1]
      jtox <- allj[j,2]
      if (pos==1) {
        jfromy <- rhist[[i]][jfromx-xlim[1]+1]
        jtoy <- rhist[[i]][jtox-xlim[1]+1]
        midy <- max(jfromy, jtoy) * 1.5
        inc <- (ylim[2] - ylim[1]) / 100 * 5
      } else {
        jfromy <- 0
        jtoy <- 0
        inc <- - (ylim[2] - ylim[1]) / 100 * 5
        midy <- inc*6
      }
      spl <- xspline(c(jfromx, jfromx, (jfromx+jtox)/2, jtox, jtox),
                     c(jfromy+inc, jfromy+2*inc, midy, jtoy+2*inc, jtoy+inc),
                     shape=c(0,1,1,1,0), draw=FALSE)
      lines(spl$x, spl$y, col=sampleColors[i], xpd=NA, lwd=.5)
      wid <- strwidth(jc[j])
      hei <- strheight(jc[j])
      tpos <- if (pos==1) max(spl$y) else min(spl$y)
      rect((jfromx+jtox)/2-wid, tpos-hei,
           (jfromx+jtox)/2+wid, tpos+hei, col="white",
           border=NA, xpd=NA)
      text((jfromx+jtox)/2, tpos, jc[j], xpd=NA)
    }
    segments(xlim[1], 0, xlim[2], 0, lty=1, lwd=.1, col=sampleColors[i])
    text(xlim[1], ylim[2], adj=c(0,0), names(reads)[i],
         xpd=NA, col=sampleColors[i])
  }

  ## Axis
  par(mar=c(0,5,0,1)+.1)
  plot(NA, type="n", xlim=xlim, ylim=ylim, ylab="", axes=FALSE,
       xlab="Genomic coordinate", xpd=NA)
  axis(1, cex.axis=.8, lwd=.4)
  print("Axis:"); print(par("fin"))
  
  ## Plot the isoforms
  isoformColorsTrans <- paste(isoformColors, sep="", "66")
  ylimgene <- c(0.5, noIso(gene)+.5)
  par(mar=c(1,5,5,1)+.1)
  plot(NA, type="n", xlim=xlim, ylim=ylimgene, xlab="", ylab="", axes=FALSE)
  print("Gene structure:") ; print(par("fin"))
  for (i in seq_along(start)) {
    segments(xlim[1], length(start)-i+1, xlim[2], length(start)-i+1,
             lty=1, lwd=.4)
    ax <- seq(xlim[1], xlim[2], length.out=40)[2:39]
    for (j in seq_along(ax)) {
      tinc <- (xlim[2]-xlim[1])/100
      xspline(c(ax[j]-tinc, ax[j], ax[j]-tinc),
              (length(start)-i+1) + c(-.1, 0, .1), shape=c(0,1,0), lwd=.4)
    }
    rect(xleft=start[[i]], xright=end[[i]], ybottom=length(start)-i+1-.25,
         ytop=length(start)-i+1+.25, col=isoformColors[i], border=NA)
    text(xlim[1], length(start)-i+1+.35, adj=c(0,0), getIso(gene)[[1]][i],
         xpd=NA, col=isoformColors[i])
  }

  if (!is.null(misoResult)) {
    breaks <- seq(0, 1, length=41)
    for (i in seq_along(misoResult)) {
      hi <- lapply(1:noIso(gene), function(j) {
        hist(misoResult[[i]]$samples[j,], breaks=breaks, plot=FALSE)
      })
      ylim <- c(0, max(sapply(hi, "[[", "counts")) * 1.3)
      par(mar=c(0,1,1,5)+.1)
      plot(NA, type="n", xlim=0:1, ylim=ylim, axes=FALSE)
      print("MISO:") ; print(par("fin"))
      sapply(seq_along(hi), function(h) {
        rect(hi[[h]]$mids-.01, 0, hi[[h]]$mids+.01, hi[[h]]$counts,
             col=isoformColorsTrans[h], border=NA)
      })
      abline(v=rowMeans(misoResult[[i]]$samples), col=isoformColors)
      axis(2, las=1, cex.axis=.8, lwd=.4)
      axis(1, at=pretty(0:1,n=4), labels=rep("", length(pretty(0:1,n=4))),
           lwd=.4)
      tl <- (par("usr")[4]-par("usr")[3])/20
      ## TODO: triangles
      segments(x0=postMean(misoResult[[i]]), y0=par("usr")[3],
               y1=par("usr")[3]-tl, col=isoformColorsTrans, lwd=3, xpd=NA,
               lend=2)
      annt <- apply(confint(misoResult[[i]]), 2,
                    function(x) { paste(round(x,2), collapse="-")})
      text(1, ylim[2], adj=c(0,1), xpd=NA, cex=.8,
           paste(round(rowMeans(misoResult[[i]]$samples), 2),
                 " [", annt, "]", sep="", collapse="\n"), col=1)
    }
    par(mar=c(0,1,0,5)+.1)
    plot(NA, type="n", xlim=0:1, ylim=0:1, xlab="", ylab="", axes=FALSE)
    title(xlab="MISO estimates", xpd=NA)
    axis(1, cex.axis=0.8, lwd=.4)
  }
  
}
