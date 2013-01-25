trainPMDHMM <-
function(m, chr.sel, nCGbin, num.cores, plot.distr=TRUE, pdfFilename){

  message("training PMD-HMM on chromosome ", chr.sel)

  indx <- as.character(seqnames(m))==chr.sel
  T <- as.numeric(values(m[indx])[, 1])
  M <- as.numeric(values(m[indx])[, 2])
  score <- calculateAlphaDistr(M, T, nCGbin, num.cores)

  
 # use parameters obtained from training on human IMR90 methylome as starting values
  J=2;
  init0 <- c(0, 1);
  P0 <- t(matrix(c(0.998297563, 0.001702437, 0.002393931, 0.997606069), nrow=J, ncol=J));
  b0 <- list(mu=c(0.3867895, 1.1690474), sigma=c(0.01649962, 0.14378640))
  startval <- hmmspec(init=init0, trans=P0, parms.emission=b0, dens.emission=dnorm.hsmm);
# train
  train <- list(x=score, N=length(score));
  startval <- hmmfit(train, startval, mstep=mstep.norm)$model

  if(plot.distr){
    x=seq(0, 3, by=0.01)
    if(!is.null(pdfFilename)){
    pdf(pdfFilename, width=5, height=5)}
    hist(score, probability=TRUE, breaks=30, xlab=sprintf("posterior mean of alpha (%s)", chr.sel), main="");lines(x, dnorm(x, mean=startval$parms.emission$mu[1], sd=sqrt(startval$parms.emission$sigma[1])), type='l', col="red");lines(x, dnorm(x, mean=startval$parms.emission$mu[2], sd=sqrt(startval$parms.emission$sigma[2])), type='l', col="green");
    if(!is.null(pdfFilename))
    dev.off()
  }
  
  startval
  
}
