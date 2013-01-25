calculateAlphaDistr <-
function(M, T, nCGbin, num.cores){

  binSize <- 0.1 # bin size to discretize alpha values
  alphas <- seq(binSize, 3, by=binSize) # discretized alphas

  logPs <- mclapply(alphas, function(a){
    as.vector(runsum(Rle(lbeta(M+a, T-M+a) - lbeta(a, a)), k=nCGbin, endrule="constant")) 
  }, mc.cores=num.cores)
  
  logPs <- do.call("cbind", logPs)
  score <- apply(logPs, 1, function(x){m=max(x); p=exp(x-m); p=p/sum(p); sum(p*alphas)})
  score

}
