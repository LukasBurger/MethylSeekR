plotPMDSegmentation <-
function(m, segs, numRegions = 1, pdfFilename=NULL, minCover = 5){

  m=m[values(m)[, 1]>=minCover]
  col.PMD="#4DAF4A"
  chrs=unique(as.character(seqnames(m)))

  height=1
  len=5000;
  

  
 # don't plot regions from chromosomes with too few CpGs
  nCGperChr=table(seqnames(m))
  indx=which(nCGperChr[chrs]<6*len)
  if(length(indx)>0){
    chrs=chrs[-indx]
  }
    
  if(!is.null(pdfFilename))
    pdf(pdfFilename, width = 15, height = 10)
  
  if(is.null(pdfFilename))
    numRegions <- 1
  
  for(ii in 1:numRegions){
    
    chr=sample(chrs, 1)
    i=as.character(seqnames(m))==chr
    s=sample(1:(length(m[i])-6*len), 1);

    par(mfrow=c(6,1), mar=c(4, 4, 3, 2))

    for(gg in 1:6){
      
      indx=(s+(gg-1)*len):(s+(gg)*len)

      PMD.sel=subsetByOverlaps(segs[values(segs)$type=="PMD"], m[i][indx])
      
      plot(start(m[i])[indx], 100*values(m[i])[indx, 2]/values(m[i])[indx, 1], pch='*', ylim=c(-15, 100), cex=0.9, ylab="methylation", xlab=sprintf("position on %s", chr))
      if(length(PMD.sel)>0){
        rect(start(PMD.sel), -15-height/4, end(PMD.sel), -15+height/4, lwd=2, col=col.PMD, border=col.PMD)
      }
      
    }
  }
  if(!is.null(pdfFilename))
    dev.off()
}
