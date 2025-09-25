\name{segmentUMRsLMRs}
\alias{segmentUMRsLMRs}
\title{
Identify UMRs and LMRs
}

\description{This function identifies hypomethylated regions and
classifies them into UMRs and LMRs.}

\usage{segmentUMRsLMRs(m, meth.cutoff = 0.5, nCpG.cutoff = 3, PMDs = NA,
pdfFilename, num.cores = 1, myGenomeSeq, seqLengths, nCpG.smoothing = 3,
minCover = 5)}

\arguments{
  
  \item{m}{GRanges object containing the methylation data.}
  
  \item{meth.cutoff}{Cut-off on methylation for calling hypomethylated
    regions.}
  
  \item{nCpG.cutoff}{Cut-off on the minimal number of CpGs for calling
    hypomethylated regions.}
  
  \item{PMDs}{GRanges object of PMDs. Set to either the return value of
    the function segmentPMDs (see example) or to NA if there are no PMDs
    (default).}
  
  \item{pdfFilename}{Name of the pdf file in which the figure is
    saved. If no name is provided (default), the figure is printed to
    the screen.}
  
  \item{num.cores}{Number of cores used for the calculations.}
  
  \item{myGenomeSeq}{Genome sequence as BSgenome object.}
  
  \item{seqLengths}{A named vector indicating the chromosome lengths of
    the genome used.}
  
  \item{nCpG.smoothing}{The number of consecutive CpGs that the
    methylation levels are averaged over.}
  
  \item{minCover}{Only CpGs with a coverage of at least minCover reads
    will be used.}

}

\value{Returns a GRanges object containing all UMRs and LMRs with the
  following metadata values: the number of CpGs with a coverage of at
  least minCover per region (nCG.segmentation), the number of CpGs in
  the DNA sequence (nCG), the total number of reads that map to CpGs in
  the region (T), the total number of read that map to CpGs without
  conversion of the C (M), the mean methylation of the segment (pmeth),
  the median methylation of the segment (median.meth) and the type
  (UMR/LMR) of region (type). The function creates a figure showing the
  classification of regions into UMRs and LMRs based on the number of
  CpGs they contain. The figure is either printed to the screen
  (default) or saved as a pdf if a filename is provided.}


\author{
  Lukas Burger lukas.burger@fmi.ch
}

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
nCpG.cutoff=n.sel, myGenomeSeq=Hsapiens, seqLengths=sLengths)

}
