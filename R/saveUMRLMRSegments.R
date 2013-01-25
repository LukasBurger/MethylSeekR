saveUMRLMRSegments <-
function(segs, GRangesFilename = NULL, TableFilename = NULL){

  # save GRanges object
  if(!is.null(GRangesFilename))
  saveRDS(segs, GRangesFilename)

  # save a tab-delimited table
  if(!is.null(TableFilename)){
    D=data.frame(chr=as.character(seqnames(segs)), start=start(segs), end=end(segs), type=values(segs)$type, nCG.segmentation=values(segs)$nCG.segmentation, nCG.seq=values(segs)$nCG, mean.meth=100*values(segs)$pmeth, median.meth=100*values(segs)$median.meth)
    write.table(D, file=TableFilename, quote=FALSE, sep="\t", row.names=FALSE)
  }
  
}
