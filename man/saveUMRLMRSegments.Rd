\name{saveUMRLMRSegments}
\alias{saveUMRLMRSegments}
\title{
Save UMR and LMR segments
}
\description{
Save UMR and LMRs segments in rds format and/or as tab-delimited file
}
\usage{
saveUMRLMRSegments(segs, GRangesFilename = NULL, TableFilename = NULL)
}
\arguments{

  \item{segs}{GRanges object containing the UMR/LMR segmentation. Return
    value of the segmentUMRsLMRs function (see example).}
  
  \item{GRangesFilename}{Filename of the GRanges object.}
  
  \item{TableFilename}{Filename of the UMR/LMR table.}
  
}

\value{No return value. Only one filename is required.}

\author{Lukas Burger lukas.burger@fmi.ch}


\examples{

library(MethylSeekR)

# get chromosome lengths
library("BSgenome.Hsapiens.UCSC.hg38")
sLengths=seqlengths(Hsapiens)

# read methylation data
methFname <- system.file("extdata", "Lister2009_imr90_hg38_chr22.tab",
package="MethylSeekR")
meth.gr <- readMethylome(FileName=methFname, seqLengths=sLengths)

FDR.cutoff <- 5 
m.sel <- 0.5 
n.sel <- 3

#segment UMRs and LMRs, assuming no PMDs
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel,
nCpG.cutoff=n.sel, num.cores=1,
myGenomeSeq=Hsapiens, seqLengths=sLengths)

#save UMRs and LMRs
saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename="UMRsLMRs.gr.rds",
 TableFilename="UMRsLMRs.tab")




}
