calculateFDRs <-
function(m, CGIs, PMDs = NA, pdfFilename=NULL, num.cores = 1, nCpG.smoothing = 3, meth.cutoffs = seq(0.3, 0.7, by=0.1), nCpG.cutoffs = seq(1, 6, by=1), minCover = 5){

  # select CpGs with minimal coverage
  m=m[values(m)[, 1]>=minCover]
  
  message("calculating false discovery rate")
  
  ## pre-allocate parameters list
  n <- length(meth.cutoffs) * length(nCpG.cutoffs)
  parameters <- vector("list", n)

  count=1
  for(meth.cutoff in meth.cutoffs){
    for(k in nCpG.cutoffs){
      parameters[[count]] <- c(meth.cutoff, k)
      count <- count+1
    }
  }

  meth.notrand <- m
  if(class(PMDs)=="GRanges"){
    message("removing PMDs for randomization")
    meth.notrand <- subsetByOverlaps(m, PMDs[values(PMDs)$type=="notPMD"])
  }

  # calculate randomized methylome
  meth.rand <- meth.notrand

  # remove CpGs in islands

  ov <- findOverlaps(meth.rand, CGIs)
  meth.rand <- meth.rand[-unique(queryHits(ov))]
  values(meth.rand) <- values(meth.rand)[sample(length(meth.rand)), ]
    
  
  res <- mclapply(parameters, function(params){
    
    
    meth.cutoff <- params[1]
    k <- params[2]
    
    # ignore chromosome boundaries for FDR calculation
    
    mean.meth <- runmean(Rle(values(meth.notrand)[, 2]/values(meth.notrand)[, 1]), k=nCpG.smoothing, endrule="constant")
    indx <- mean.meth < meth.cutoff
    nSeg <- sum(runValue(indx)==TRUE & runLength(indx)>=k)

    mean.meth=runmean(Rle(values(meth.rand)[, 2]/values(meth.rand)[, 1]), k=nCpG.smoothing, endrule="constant")
    indx <- mean.meth < meth.cutoff
    nSeg.rand <- sum(runValue(indx)==TRUE & runLength(indx)>=k)

    c(nSeg, nSeg.rand)
    
  }, mc.cores=num.cores)
  
  FDRs=100*sapply(res, function(x){x[2]/x[1]})

  tmp=matrix(NA, nrow=length(meth.cutoffs), ncol=length(nCpG.cutoffs))
  rownames(tmp)=meth.cutoffs
  colnames(tmp)=nCpG.cutoffs

  
  count=1
  for(meth.cutoff in meth.cutoffs){
    for(k in nCpG.cutoffs){
      tmp[as.character(meth.cutoff), as.character(k)]=FDRs[count]
      count=count+1
    }
  }
  FDRs=tmp
  
  numSegments=sapply(res, function(x){x[1]})

  tmp=matrix(NA, nrow=length(meth.cutoffs), ncol=length(nCpG.cutoffs))
  rownames(tmp)=meth.cutoffs
  colnames(tmp)=nCpG.cutoffs

  count=1
  for(meth.cutoff in meth.cutoffs){
    for(k in nCpG.cutoffs){
      tmp[as.character(meth.cutoff), as.character(k)]=numSegments[count]
      count=count+1
    }
  }
  numSegments=tmp

  rownames(FDRs)=as.character(100*as.numeric(rownames(FDRs)))
  rownames(numSegments)=as.character(100*as.numeric(rownames(numSegments)))

  if(!is.null(pdfFilename)){
    pdf(pdfFilename, width = 9, height=4.5)}
  
  par(mfrow=c(1,2))
  barplot(pmin((t(FDRs)), 20), beside=TRUE, ylab="FDR (%)", ylim=c(0,20), xlab="methylation cut-off (%)")
  barplot(t(numSegments), beside=TRUE, ylab="number of segments", xlab="methylation cut-off (%)")
  legend("topleft", legend=paste(colnames(FDRs), c("CpG", rep("CpGs", ncol(FDRs)-1)), sep=" "), fill=grey.colors(ncol(FDRs)), bty="n")
  if(!is.null(pdfFilename))
    dev.off()

  rownames(FDRs)=as.character(as.numeric(rownames(FDRs))/100)
  rownames(numSegments)=as.character(as.numeric(rownames(numSegments))/100)
  list(FDRs=FDRs, numSegments=numSegments)
}
