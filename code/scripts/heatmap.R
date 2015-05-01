
library(SCM)
library(RColorBrewer)

data <- loadExprAbund(transformation="log")
allData <- cbind(data$data1, data$data2)[, sort(c(colnames(data$data1), 
                                                  colnames(data$data2)))]

## Merge ingolia
colnames(allData) <- gsub("ingolia([123])\\.", "ingolia.\\1",
                          perl = TRUE, colnames(allData))

type <- c(velc="SAGE", roth="MA", holstege="MA", causton="MA",
          dudley="cMA", miura="cPCR", nagalakshmi="RS",
          yassour="RS", lipson="RS", lipson2="MA", ingolia="RS",
          mackay="cMA", garcia="cMA", pelechano="cMA", futcher="gel",
          gygi="gel", newman="gfp", degodoy="ms", lee="ms", lu="ms",
          nagaraj="ms", peng="ms", thakur="ms", washburn="ms",
          ghaem="wblot")

experiment <- sapply(strsplit(colnames(allData), ".", fixed = TRUE), "[[", 2)
exp_type <-type[experiment]

allData <- allData[, order(exp_type)]

allData <- allData[, c(grep("^expr", colnames(allData)),
                       grep("^abund", colnames(allData))) ]

CC <- cor(allData, use="pairwise.complete.obs")

ticks <- which(!duplicated(sapply(strsplit(colnames(allData), "\\."),
                                  "[", 2)))[-1]

scaler <- function(range, domain) {
  stopifnot(length(range)==2, length(domain)==2,
            is.numeric(range), is.numeric(domain),
            range[1] < range[2], domain[1] < domain[2])
  function(x, invert=FALSE) {
    (x - domain[1]) / (domain[2]-domain[1]) * (range[2]-range[1]) + range[1]
  }
}

myheatmap <- function(mat, xbreaks=numeric(), ybreaks=numeric(),
                      xbreaksize=.1, ybreaksize=xbreaksize,
                      rectborder=par("fg"), lwd=0.01,
                      colors=brewer.pal(9, "Blues"),
                      colrange=range(mat), xlab="", ylab="",
                      key=TRUE, mar=c(5,9,5,1)+.1, asp=TRUE,
                      sep_color = "white", ...) {

  xbreaksize <- rep(xbreaksize, length=length(xbreaks))
  ybreaksize <- rep(ybreaksize, length=length(ybreaks))
  
  ## colors
  breaks <- seq(colrange[1], colrange[2], length=length(colors)+1)
  rectcol <- as.integer(cut(mat, breaks=breaks))

  ## rectange positions
  xleft <- matrix(rep(0:(ncol(mat)-1), nrow(mat)),
    nrow=nrow(mat), byrow=TRUE)
  ybottom <- matrix(rep(0:(nrow(mat)-1), each=ncol(mat)),
                    nrow=nrow(mat), byrow=TRUE)
  xright <- matrix(rep(1:ncol(mat), nrow(mat)),
    nrow=nrow(mat), byrow=TRUE)
  ytop <- matrix(rep(1:nrow(mat), each=ncol(mat)),
                 nrow=nrow(mat), byrow=TRUE)

  ## add breaks
  shiftx <- function(x, from, by) {
    x[,from:ncol(x)] <- x[, from:ncol(x)] + by
    x
  }
  shifty <- function(x, from, by) {
    x[from:nrow(x),] <- x[from:nrow(x),] + by
    x
  }
  for (i in seq_along(xbreaks)) {
    xleft <- shiftx(xleft, from=xbreaks[i], by=xbreaksize[i])
    xright <- shiftx(xright, from=xbreaks[i], by=xbreaksize[i])
  }
  for (i in seq_along(ybreaks)) {
    ytop <- shifty(ytop, from=ybreaks[i], by=ybreaksize[i])
    ybottom <- shifty(ybottom, from=ybreaks[i], by=ybreaksize[i])
  }

  ## Leave out space for key
  realxlim <- xlim <- range(c(xleft, xright))
  ylim <- range(c(ytop, ybottom))
  if (key) {
    realxlim <- c(xlim[1], xlim[2]+5)
  }
  
  ## empty plot
  par(mar = mar)
  plot(NA, type="n", xlim=realxlim, ylim=ylim,
       frame=FALSE, axes=FALSE, xlab="", ylab="", asp=asp, ...)

  ## plot rectangles
  rect(xleft=xleft, ybottom=max(ytop) - ybottom,
       xright=xright, ytop=max(ytop) - ytop, 
       col=colors[rectcol], lwd=lwd, border=rectborder)

  ## make breaks grey
  greyleft <- xright[1,][-nrow(xright)]
  greyright <- xleft[1,][-1]
  greyplot <- which(greyleft != greyright)
  rect(xleft=greyleft[greyplot], xright=greyright[greyplot],
       ybottom=ylim[1], ytop=ylim[2], col=sep_color, border=NA)

  greytop <- ytop[,1][-nrow(ytop)]
  greybottom <- ybottom[,1][-1]
  greyplot2 <- which(greytop != greybottom)
  rect(xleft=xleft[1,1], xright=xright[1,ncol(xright)],
       ybottom=max(ytop)-greybottom[greyplot2],
       ytop=max(ytop)-greytop[greyplot2], col=sep_color, border=NA)
  
  ## color key
  if (key) {
    kxleft <- realxlim[2]-2
    kxright <- realxlim[2]
    kybreaks <- seq(ylim[1], ylim[2], length=length(colors)+1)
    kytop <- kybreaks[-1]
    kybottom <- kybreaks[-length(kybreaks)]
    rect(xleft=kxleft, xright=kxright, ytop=kytop, ybottom=kybottom,
         border=colors, col=colors)
    rect(xleft=min(kxleft), xright=max(kxright),
         ytop=max(kytop), ybottom=min(kybottom), lwd=.1)
    keyscaler <- scaler(range=ylim, domain=colrange)
    labs <- pretty(breaks)
    text(realxlim[2], keyscaler(labs)[-c(1,length(labs))],
         labs[-c(1,length(labs))], pos=4, xpd=NA)
    text(realxlim[2], keyscaler(labs)[1], labs[1], pos=4, adj=c(0,1), xpd=NA)
    text(realxlim[2], keyscaler(labs)[length(labs)],
         labs[length(labs)], pos=4, adj=c(0,1), xpd=NA)
  }

  ## return some info
  list(xlim=xlim, ylim=ylim, xleft=xleft, xright=xright, ytop=ytop,
       ybottom=ybottom)
}

pdf("cor-heat.pdf", width=10, height=8)
first.abund <- grep("^abund", colnames(CC))[1]
breaks <- c(ticks, first.abund)
breaksize <- c(rep(.1, length(ticks)), .5)
heat <- myheatmap(CC, colrange=c(0,1), rectborder=NA,
                  xbreaks=breaks, ybreaks=breaks,
                  xbreaksize=breaksize, ybreaksize=breaksize,
                  mar=c(3,2,1,1))
rect(xleft=heat$xleft[1,1], xright=heat$xright[1, first.abund-1],
     ytop=heat$ytop[nrow(CC),1]+1, ybottom=heat$ytop[nrow(CC),1]+.3, 
     col="darkgrey", border="darkgrey", xpd=NA)
rect(xleft=heat$xleft[1,first.abund], xright=heat$xright[1,ncol(CC)],
     ytop=heat$ytop[nrow(CC),1]+1, ybottom=heat$ytop[nrow(CC),1]+.3,
     col="lightgrey", border="lightgrey", xpd=NA)
text((heat$xleft[1,1] + heat$xright[1,first.abund-1]) / 2,
     heat$ytop[nrow(CC),1]+2, adj=c(1/2,1/2), "Gene expression", xpd=NA)
text((heat$xleft[1,first.abund] + heat$xright[1,ncol(CC)]) / 2,
     heat$ytop[nrow(CC),1]+2, adj=c(1/2,1/2), "Protein abundance", xpd=NA)

cap <- function(s) {
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep="")
}

pn0 <- sub("\\.[^\\.]*$", "", sub("^[^\\.]*\\.", "", colnames(CC)))
pn <- cap(pn0)
pn[duplicated(pn)] <- ""
pn <- sub("([0-9])", " #\\1", pn)

rcol <- rep("lightgrey", length(pn))
pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
         "#000000")
rcol[type[pn0]=="SAGE"] <- pal[1]
rcol[type[pn0]=="MA"]   <- pal[2]
rcol[type[pn0]=="cMA"]  <- pal[3]
rcol[type[pn0]=="RS"]   <- pal[4]
rcol[type[pn0]=="cPCR"] <- pal[5]
rcol[type[pn0]=="gel"] <- pal[6]
rcol[type[pn0]=="gfp"] <- pal[7]
rcol[type[pn0]=="ms"] <- pal[8]
rcol[type[pn0]=="wblot"] <- pal[9]

rect(xleft=heat$xleft[1,1]-1, xright=heat$xleft[1,1]-.3,
     ytop=max(heat$ytop) - heat$ytop[,1] +
     ifelse(c(pn[-1], "foo") =="", 0, .1),
     ybottom=max(heat$ytop) - heat$ybottom[,1] - ifelse(pn=="", 0, .1),
     border=rcol, col=rcol)

text(heat$xleft[1,1]-.5,
     max(heat$ytop) - (heat$ytop[1:nrow(CC)]+heat$ybottom[1:nrow(CC)])/2 -
     ifelse(c(pn[-1], "foo")=="", .5, 0),
     pos=2, pn, cex=.6, xpd=NA)

legend(0, -1, fill=pal, border=pal,
       c("SAGE", "Comm. microarray", "Custom microarray", "RNA-Seq", "cPCR", "2d gel",
         "GFP", "LC MS/MS", "Western blot"), ncol=5, bty = 'n', cex = .6, xpd = NA,
       pt.cex = 1.2)

dev.off()
