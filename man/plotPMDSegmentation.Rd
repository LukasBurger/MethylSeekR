\name{plotPMDSegmentation}
\alias{plotPMDSegmentation}
\title{
Plotting the PMD Segmentation 
}
\description{
This function generates a figure showing the PMD segmentation in a randomly chosen region.
}
\usage{
plotPMDSegmentation(m, segs, numRegions = 1, pdfFilename=NULL, minCover = 5)
}
\arguments{
  \item{m}{
    GRanges object containing the methylation data.
  }
  \item{segs}{
    GRanges object containing the PMD segmentation. Return value of the
    segmentPMDs function (see example).
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
}

\value{

  No return value. The function creates a figure showing the inferred
  segmentation for a randomly chosen region. The figure is either
  printed to the screen (default) or saved as a pdf if a filename is
  provided. If a filename (pdfFilename) is provided, several regions
  (set via the numRegions argument) can be plotted and saved in a
  multi-page pdf file. The randomly chosen region that is displayed is
  broken up into 6 panels and in each panel, the raw (ie unsmoothed)
  methylation levels of all CpGs with a minimal coverage of 5 reads are
  shown. PMDs are indicated as green bars, extending over the entire
  PMD.
  
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

#segment PMDs
PMDsegments.gr <- segmentPMDs(m=meth.gr, chr.sel="chr22", seqLengths=sLengths)

#plot PMD segmentation examples
plotPMDSegmentation(m=meth.gr, segs=PMDsegments.gr, numRegions=1)


}
