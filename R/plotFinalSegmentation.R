plotFinalSegmentation <-
function(m, segs, PMDs = NA, meth.cutoff, numRegions = 1, pdfFilename = NULL,  minCover = 5, nCpG.smoothing = 3){

  m=m[values(m)[, 1]>=minCover]
  
  cols.seg=c("#377EB8", "#E41A1C")
  names(cols.seg)=c("UMR", "LMR")

  col.PMD="#4DAF4A"
  
  pch.seg <- c(LMR=24, UMR=22)
  chrs=unique(as.character(seqnames(m)))

  height=1
  len=5000;

  # don't plot regions from chromosomes with too few CpGs
  nCGperChr=table(seqnames(m))
  indx=which(nCGperChr[chrs]<3*len)
  if(length(indx)>0){
    chrs=chrs[-indx]
  }
  

  exist.PMDs <- FALSE
  if(class(PMDs)=="GRanges"){
    exist.PMDs <- TRUE # are there any PMDs ?
   }

  if(!is.null(pdfFilename))
    pdf(pdfFilename, width = 15, height = 8)

  if(is.null(pdfFilename))
    numRegions <- 1
  
  for(ii in 1:numRegions){

    chr=sample(chrs, 1)
    i=seqnames(m)==chr
    s=sample(1:(length(m[i])-3*len), 1);

    
    par(mfrow=c(6,1))

    for(gg in 1:3){
      
      indx=(s+(gg-1)*len):(s+(gg)*len)
      
      seg.sel=subsetByOverlaps(segs, m[i][indx])
      if(length(seg.sel)>0){
        mysegcol <- cols.seg[as.character(values(seg.sel)$type)]
        mysegpch <- pch.seg[as.character(values(seg.sel)$type)]
      }

      PMD.sel=NA
      if(exist.PMDs){
        PMD.sel=subsetByOverlaps(PMDs[values(PMDs)$type=="PMD"], m[i][indx])
      }


      add.mar <- 2
      
      par(mar=c(3-add.mar, 4, 1+add.mar, 2))
      plot(start(m[i])[indx], 100*values(m[i])[indx, 2]/values(m[i])[indx, 1], pch='*', ylim=c(-18, 100), cex=0.9, ylab="methylation", axes=FALSE, xlab="")
      axis(2, at=c(0, 50, 100))
      box()
      text(par("usr")[1]+par("cxy")[1]*0.2, par("usr")[3]+par("cxy")[2]*0.3, adj=c(0, 0), label="raw") 
      
      if(exist.PMDs & length(PMD.sel)>0){
        rect(start(PMD.sel), -18-height/4, end(PMD.sel), -18+height/4, lwd=2, col=col.PMD, border=col.PMD)
      }
      if(length(seg.sel)>0){
        points(x=mid(ranges(seg.sel)), y=rep(-18,length(seg.sel)), pch=mysegpch, cex=1.3, col="black", lwd=0.5, bg=mysegcol)
      }

      par(mar=c(2+add.mar, 4, 2-add.mar, 2))
      plot(start(m[i])[indx],100*as.vector(runmean(Rle(values(m[i])[indx, 2]/values(m[i])[indx, 1]), nCpG.smoothing, endrule="constant")), pch='*', ylim=c(-18, 100), cex=0.9, ylab="methylation", xlab=sprintf("position on %s", chr), axes=FALSE);
      axis(1)
      axis(2, at=c(0,50,100))
      box()
      text(par("usr")[1]+par("cxy")[1]*0.2, par("usr")[3]+par("cxy")[2]*0.3, adj=c(0, 0), label="smoothed") 
      
      if(exist.PMDs & length(PMD.sel)>0){
        rect(start(PMD.sel), -18-height/4, end(PMD.sel), -18+height/4, lwd=2, col=col.PMD, border=col.PMD)
      }
      if(length(seg.sel)>0){
        points(x=mid(ranges(seg.sel)), y=rep(-18,length(seg.sel)), pch=mysegpch, cex=1.3, col="black", lwd=0.5, bg=mysegcol)
      }
      abline(h=100*meth.cutoff, lty=5, col="darkgrey")


      
    }
  }
  if(!is.null(pdfFilename))
    dev.off()
 
}
