PMDviterbiSegmentation <-
function(m, hmm.model, nCGbin, num.cores){

  message("performing viterbi segmentation")

  chrs=unique(as.character(seqnames(m)))
  
  y.list=mclapply(chrs, function(chr.sel){

    indx=as.character(seqnames(m))==chr.sel;
    T <- as.numeric(values(m[indx])[, 1])
    M <- as.numeric(values(m[indx])[, 2])
    score <- calculateAlphaDistr(M, T, nCGbin, num.cores)
    train=list(x=score, N=length(score));
    y=predict(hmm.model, train);

    #remove regions that are too short
    ttt=Rle(y$s)
    min.len=nCGbin

    # first take regions that are PMD, but too short and make them nonPMD
    indx=runLength(ttt)<=min.len & runValue(ttt)==2;
    runValue(ttt)[indx]=1;
    # now vice versa
    indx=runLength(ttt)<=min.len & runValue(ttt)==1;
    runValue(ttt)[indx]=2;
    y$s=as.vector(ttt);
    y$score=score
    
    return(y)

  }, mc.cores=num.cores);

  names(y.list)=chrs
  y.list

}
