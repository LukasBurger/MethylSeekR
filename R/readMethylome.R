readMethylome <-
function(FileName, seqLengths, format="text"){

  message("reading methylome data")
  if(format=="GRanges"){
    meth <- readRDS(FileName)
  }
  else if(format=="text"){
    dat <- scan(FileName, what=list(chr=character(0), pos=integer(0), T=double(0), M=double(0)), sep="\t")
    meth <- GRanges(seqnames=dat$chr, ranges=IRanges(start=dat$pos, width=1), strand='*', T=dat$T, M=dat$M, seqlengths=seqLengths)
  }
  else{
    stop("unknown format")
  }

  mean.cov=mean(values(meth)[, 1])
  if(mean.cov < 10)
  warning(sprintf("We do not recommend the use of MethylSeekR\nfor methylomes with mean coverage < 10X\n(mean coverage of CpGs with at least one read: %.1f)", mean.cov))
  # order by chromosome and CpG position
  meth <- meth[order(as.vector(seqnames(meth)), start(meth))]
  meth
  
}
