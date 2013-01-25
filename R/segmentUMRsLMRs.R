segmentUMRsLMRs <-
function(m, meth.cutoff = 0.5, nCpG.cutoff = 3, PMDs = NA, pdfFilename = NULL, num.cores = 1, myGenomeSeq, seqLengths, nCpG.smoothing = 3, minCover = 5){


  nCG.classification <- 30 # cut-off for the separation of LMRs and UMRs
  
  message("identifying UMRs and LMRs")
 # select CpGs with minimal coverage
  m=m[values(m)[, 1]>=minCover]

  chrs=unique(as.character(seqnames(m)))
    
  res <- mclapply(chrs, function(chr){
    
    sel <- which(as.character(seqnames(m))==chr)
    mean.meth <- runmean(Rle(values(m)[sel, 2]/values(m)[sel, 1]), k=nCpG.smoothing, endrule="constant")

    indx <- mean.meth < meth.cutoff
    
    # remove segments that are too short
    # first hypomethylated regions
    runValue(indx)[runLength(indx)<nCpG.cutoff & runValue(indx)==TRUE]=FALSE
    # then methylated regions
    runValue(indx)[runLength(indx)<nCpG.cutoff & runValue(indx)==FALSE]=TRUE

    tmp.ids <- rep(1:length(runLength(indx)), runLength(indx))
    tmp <- split(1:length(sel), tmp.ids)

      
    tmp <- tmp[runValue(indx)==TRUE]

    if(length(tmp)>0){

      coords <- cbind(sapply(tmp, min), sapply(tmp, max))
      starts <- round((start(m)[sel[pmax(1, coords[, 1]-1)]]+start(m)[sel[coords[, 1]]])/2) # take middle between two CpGs (if at the beginning or end of chromosome, don't extend
      ends <- round((start(m)[sel[coords[, 2]]]+start(m)[sel[pmin(length(sel), coords[, 2]+1)]])/2)
          
      hmr.gr=GRanges(seqnames=unique(seqnames(m[sel])), strand="*", ranges=IRanges(starts, ends), seqlengths=seqLengths)

    }
    else{
      
      hmr.gr=GRanges(, seqlengths=seqLengths)

    }
    hmr.gr
  }, mc.cores=num.cores)

  segments.gr=do.call(c, unname(res))

  # remove segments that lie in PMDs if there are any
  if(class(PMDs)=="GRanges"){
    segments.gr=subsetByOverlaps(segments.gr, PMDs[values(PMDs)$type=="notPMD"])
  }
  
  nCG=vcountPattern("CG", getSeq(myGenomeSeq, resize(segments.gr, width(segments.gr), fix="start"), as.character=FALSE))#
  ov <- findOverlaps(m, segments.gr)
  T=tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), sum)
  M=tapply(values(m[queryHits(ov)])[, 2], subjectHits(ov), sum)
  nCG.segmentation=tapply(values(m[queryHits(ov)])[, 1], subjectHits(ov), length)
  median.meth=tapply(as.vector(runmean(Rle(values(m[queryHits(ov)])[, 2]/values(m[queryHits(ov)])[, 1]), nCpG.smoothing, endrule="constant")), subjectHits(ov), median)
  median.meth=pmax(0, median.meth) # conversion to Rle causes rounding errors that sometimes give very small negative numbers

  if(!all.equal(as.numeric(names(T)), 1:length(segments.gr))){
      message("error in calculating methylation levels for PMDs")
    }
        
  
  type=c("UMR", "LMR")[as.integer(nCG<nCG.classification)+1]
  values(segments.gr)=DataFrame(nCG.segmentation, nCG, T, M, pmeth=M/T, median.meth=median.meth, type)

  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  if(!is.null(pdfFilename)){
    pdf(pdfFilename, height=5, width=5)}
  smoothScatter(log2(values(segments.gr)$nCG), 100*values(segments.gr)$median.meth, colramp=jet.colors, xlab="log2 number of CpGs in segment", ylab="median methylation (%)");abline(v=log2(nCG.classification), lty=5)
  if(!is.null(pdfFilename))
    dev.off()
  
  segments.gr
    
}
