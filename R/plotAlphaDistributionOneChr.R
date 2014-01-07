plotAlphaDistributionOneChr <-
function(m, chr.sel, pdfFilename = NULL, num.cores = 1, nCGbin = 101){

  message("determining alpha distribution for chromosome: ", chr.sel)
  
  indx <- as.character(seqnames(m))==chr.sel
  if(sum(indx)<nCGbin)
    stop(sprintf("Error: less than %d covered CpGs on chromosome %s", nCGbin, chr.sel))
  
  T <- as.numeric(values(m[indx])[, 1])
  M <- as.numeric(values(m[indx])[, 2])

  score <- calculateAlphaDistr(M, T, nCGbin, num.cores)

  if(!is.null(pdfFilename)){
    pdf(pdfFilename, width=5, height=5)
  }
  
  hist(score, probability=TRUE, breaks=30, xlab=sprintf("posterior mean of alpha (%s)", chr.sel), main="")
  if(!is.null(pdfFilename))
  dev.off()

}
