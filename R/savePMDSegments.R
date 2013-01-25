savePMDSegments <-
function(PMDs, GRangesFilename = NULL, TableFilename = NULL){

  # save as GRanges object
  if(!is.null(GRangesFilename))
  saveRDS(PMDs, GRangesFilename)

  # save as tab-delimited table

  if(!is.null(TableFilename)){
    indx=values(PMDs)$type=="PMD"
    D=data.frame(chr=as.character(seqnames(PMDs))[indx], start=start(PMDs)[indx], end=end(PMDs)[indx])
    write.table(D, file=TableFilename, quote=FALSE, sep="\t", row.names=FALSE)
  }
  
}
