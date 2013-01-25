createGRangesObjectPMDSegmentation <-
function(m, y.list, num.cores, seqLengths){

  message("creating GRanges object")
  chrs=unique(as.character(seqnames(m)))

  segList <- mclapply(chrs, function(chr.sel) {

    indx=as.character(seqnames(m))==chr.sel
    n <- sum(indx)
    mids <- round(0.5*(start(m[indx])[-length(m[indx])] + start(m[indx])[-1]))
    segCG <- Rle(y.list[[chr.sel]]$s)
    segNt <- Rle(lengths=c(diff(c(1,mids)),seqLengths[chr.sel]-mids[n-1]+1), values=y.list[[chr.sel]]$s)
    segChr <- GRanges(seqnames=chr.sel, IRanges(start=c(1,cumsum(runLength(segNt))[-nrun(segNt)]+1), end=cumsum(runLength(segNt))), strand="*", type=c("notPMD","PMD")[runValue(segNt)], nCG=runLength(segCG), seqlengths=seqLengths[chrs])
    segChr

  }, mc.cores=num.cores);

  segments <- do.call(c, unname(segList))
  segments

}
