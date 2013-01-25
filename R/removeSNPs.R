removeSNPs <-
function(m, snps){

  message("removing SNPs")
  width(m) <- 2
  ov <- findOverlaps(m, snps)
  m <- m[-unique(queryHits(ov))]
  n=length(unique(queryHits(ov)))
  message(sprintf("%d (%.1f%s) CpGs removed", n, 100*(n/length(m)), "%"))
  width(m) <- 1
  m

}
