segmentPMDs <-
function(m, chr.sel, pdfFilename = NULL, seqLengths, num.cores = 1, nCGbin = 101){ 

  hmm.model <- trainPMDHMM(m, chr.sel, nCGbin, num.cores, plot.distr=TRUE, pdfFilename)
  y.list <- PMDviterbiSegmentation(m, hmm.model, nCGbin, num.cores)
  segments <- createGRangesObjectPMDSegmentation(m, y.list, num.cores, seqLengths)
  segments
   
 }
