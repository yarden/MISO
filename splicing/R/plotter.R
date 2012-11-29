
plotIso <- function(geneStructure, gene=1, xlab="", ylab="",
                    col=rainbow_hcl(noIso(geneStructure)[gene]),
                    mar=c(0,0,2,0)+.1, axes=FALSE, ...) {

  require(colorspace)
  
  noIso <- noIso(geneStructure)[gene]  
  col <- rep(col, length.out=noIso)

  start <- getExonStart(geneStructure, gene)
  end <- getExonEnd(geneStructure, gene)
  xlim <- range(unlist(start), unlist(end))
  
  ylimgene <- c(0.5, noIso+.5)
  par(mar=mar)
  plot(NA, type="n", xlim=xlim, ylim=ylimgene, xlab=xlab,
       ylab=ylab, axes=axes, ...)

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
         ytop=length(start)-i+1+.25, col=col[i], border=NA)
    text(xlim[1], length(start)-i+1+.35, adj=c(0,0),
         getIso(geneStructure)[[gene]][i],
         xpd=NA, col=col[i])
  }
}

plotIsoSize <- function(geneStructure, gene=1,
                        geneBaseHeight=0.5, geneIsoHeight=0.5) {
  c(width=5,
    height=noIso(geneStructure)[gene] * geneIsoHeight + geneBaseHeight)
}

plotIsoPDF <- function(geneStructure, file, gene=1, mar=c(0,0,2,0), ...) {

  pdf(tmp <- tempfile())
  plot.new()
  size <- plotIsoSize(geneStructure, gene)
  dev.off()
  unlink(tmp)
  
  pdf(file, width=size[1], height=size[2])
  par(mar=mar)
  plotIso(geneStructure, gene, ...)
  dev.off()
}

plotMISO <- function(misoResult,
                     col=rainbow_hcl(noIso(misoResult$geneStructure)),
                     legend=c("topright", "topleft", "rightmargin", "none"),
                     xlab="Isoform ratio", ylab="Relative frequency",
                     frame=FALSE, axes=TRUE, cex.axis=0.8, lwd.axis=0.4,
                     las=1, ...) {

  require(colorspace)

  legend <- match.arg(legend)
  
  col <- rep(col, length.out=noIso(misoResult$geneStructure))
  colTrans <- paste(col, sep="", "66")
  
  breaks <- seq(0, 1, length=41)

  hi <- lapply(1:noIso(gene), function(j) {
    hist(misoResult$samples[j,], breaks=breaks, plot=FALSE)
  })
  ylim <- c(0, max(sapply(hi, function(x) x$counts/sum(x$counts))) * 1.3)
  plot(NA, type="n", xlim=0:1, ylim=ylim, frame=frame, axes=FALSE,
       xlab=xlab, ylab=ylab, ...)
  if (axes) {
    axis(2, las=las, cex.axis=cex.axis, lwd=lwd.axis)
    axis(1, las=las, lwd=lwd.axis, cex.axis=cex.axis)
  }
  sapply(seq_along(hi), function(h) {
    rect(hi[[h]]$mids-.01, 0, hi[[h]]$mids+.01,
         hi[[h]]$counts/sum(hi[[h]]$counts),
         col=colTrans[h], border=NA)
  })
  abline(v=rowMeans(misoResult$samples), col=col)
  tl <- (par("usr")[4]-par("usr")[3])/20
  ## TODO: triangles
  segments(x0=postMean(misoResult), y0=par("usr")[3],
           y1=par("usr")[3]-tl, col=colTrans, lwd=3, xpd=NA,
           lend=2)

  if (legend != "none") {
    annt <- apply(confint(misoResult), 2,
                  function(x) { paste(round(x,2), collapse="-")})
    if (legend == "topright") {
      xx <- 1
      yy <- ylim[2]-(1:noIso(gene)-1)*strheight("[")*1.2
      adj <- c(1,1)
    } else if (legend == "topleft") {
      xx <- 0
      yy <- ylim[2]-(1:noIso(gene)-1)*strheight("[")*1.2
      adj <- c(0,1)
    } else if (legend == "rightmargin") {
      xx <- 1
      yy <- ylim[2]-(1:noIso(gene)-1)*strheight("[")*1.2
      adj <- c(0,1)
    }
    text(xx, yy, adj=adj, xpd=NA, cex=.8,
         paste(round(rowMeans(misoResult$samples), 2),
               " [", annt, "]", sep=""), col=col)
  }
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
  plotIso(gene, mar=c(1,5,5,1)+.1, col=isoformColors)
  
  if (!is.null(misoResult)) {
    lapply(misoResult, function(ms) {
      par(mar=c(0,1,1,5)+.1)
      plotMISO(ms, col=isoformColors, legend="rightmargin", axes=FALSE)
      axis(2, lwd=0.4, cex.axis=0.8, las=1)
      axis(1, lwd=0.4, cex.axis=0.8, las=1,
           at=pretty(0:1), labels=rep("", length(pretty(0:1))))
    })
    par(mar=c(0,1,0,5)+.1)
    plot(NA, type="n", xlim=0:1, ylim=0:1, xlab="", ylab="", axes=FALSE)
    title(xlab="MISO estimates", xpd=NA)
    axis(1, cex.axis=0.8, lwd=.4)
  }  
}
