readSNPTable <-
function(FileName, seqLengths, format="text"){

  message("reading SNP table")
  
  if(format=="GRanges"){
    snp <- readRDS(FileName)
  }
  else if(format=="text"){
    dat <- scan(FileName, what=list(chr=character(0), pos=integer(0)), sep="\t")
    snp <- GRanges(seqnames=dat$chr, ranges=IRanges(start=dat$pos, width=1), strand='*', seqlengths=seqLengths)
  }
  else{
    stop("unknown format")
  }

  # order by chromosome and CpG position
  snp <- snp[order(as.vector(seqnames(snp)), start(snp))]
  snp
  
}
