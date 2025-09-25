\name{plotFinalSegmentation}
\alias{plotFinalSegmentation}
\title{
Plotting final segmentation
}
\description{
This function plots the final segmentations, including PMDs, UMRs
and LMRs.
}
\usage{
plotFinalSegmentation(m, segs, PMDs = NA, meth.cutoff, numRegions = 1,
pdfFilename = NULL, minCover = 5, nCpG.smoothing = 3)
}
\arguments{
  \item{m}{
    GRanges object containing the methylation data.
    }
  \item{segs}{
    GRanges object containing the UMR/LMR segmentation. Return value of the
    segmentUMRsLMRs function (see example).
  }
  \item{PMDs}{
    GRanges object of PMDs. Set to either the return value of the
    function segmentPMDs (see example) or to NA if there are no PMDs.
  }
  \item{meth.cutoff}{
    Cut-off on methylation for calling hypomethylated regions.
  }
  \item{numRegions}{
    The number of randomly chosen regions to be plotted. The default (1)
    can only be changed if a pdfFilename is provided (see below).
  }
  \item{pdfFilename}{
    Name of the pdf file in which the figure is saved. If no name is
    provided (default), the figure is printed to the screen.
  }
  \item{minCover}{
    Only CpGs with a coverage of at least minCover reads will be used
    for plotting.
  }
  \item{nCpG.smoothing}{
    The number of consecutive CpGs that the methylation levels are
    averaged over.
  }
}
\value{

  No return value. The function creates a figure showing the inferred
  segmentation for a randomly chosen region. The figure is either
  printed to the screen (default) or saved as a pdf if a filename is
  provided. If a filename (pdfFilename) is provided, several regions
  (set via the numRegions argument) can be plotted and saved in a
  multi-page pdf file. The randomly chosen region that is displayed is
  broken up into 3 pairs of panels, where in each pair the same region
  is shown twice, once with raw methylation levels (top) and once with
  methylation levels smoothed over 3 consecutive CpGs (bottom). In both
  cases only CpGs with a coverage of at least minCover reads are
  shown. The raw data best illustrates the disordered nature of
  methylation levels in PMDs, whereas the smoothed methylation levels
  more clearly show UMRs and LMRs. In all figures, UMRs are shown as
  blue squares (placed at the middle of the identified segment), LMRs as
  red triangles (placed at the middle of the identified segment) and
  PMDs as green bars (extending over the entire PMD). The cut-off on
  methylation (meth.cutoff) to determine UMRs and LMRs is shown as a
  grey dashed line.

  

}
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

#plot final segmentation, assuming no PMDs
plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr, meth.cutoff=m.sel, numRegions=1)


}
