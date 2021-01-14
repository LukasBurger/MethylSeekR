\name{calculateFDRs}
\alias{calculateFDRs}
\title{
Determine false discovery rate 
}
\description{
This function calculates the false discovery rate (FDR) based on randomized
Bis-seq data.
}
\usage{
calculateFDRs(m, CGIs, PMDs = NA, pdfFilename=NULL, num.cores = 1,
nCpG.smoothing = 3, meth.cutoffs = seq(0.3, 0.7, by=0.1), nCpG.cutoffs =
seq(1, 6, by=1), minCover = 5)
}
\arguments{
  \item{m}{
    GRanges object containing the methylation data.
  }
  \item{CGIs}{
    A GRanges object of CpG island coordinates. All CpGs overlapping CpG islands
    will be removed for the randomization.
  }
  \item{PMDs}{
    The GRanges object of the PMDs. Set to either the return value of the
    function segmentPMDs (see example) or to NA (default) if there are no PMDs.
  }
  \item{pdfFilename}{
    Name of the pdf file in which the figure is saved. If no name is
    provided (default), the figure is printed to the screen.
  }
  \item{num.cores}{
    Number of cores used for the calculations.
  }
  \item{nCpG.smoothing}{
   The number of consecutive CpGs that the methylation levels are
   averaged over.  
  }
  \item{meth.cutoffs}{
   A vector containing the cut-offs in methylation for which the FDR
   should be calculated. Numbers must be between 0 and 1. 
  }
  \item{nCpG.cutoffs}{
    A vector containing the cut-offs on the minimal number of CpGs per
    region for which the FDR should be calculated.
  }
  \item{minCover}{
    Only CpGs with a coverage of at least minCover reads will be used.
  }
}
\value{
  A list containing a matrix with FDR values and a matrix with the
  number of inferred segments for each methylation cut-off (rows)
  and each cut-off on the minimal number of CpGs per region
  (columns). The function creates a figure showing the relationship
  between the methylation cut-off, the cut-off on the minimal number of
  CpGs per region, the number of inferred segments and the FDR. The
  figure is either printed to the screen (default) or saved as a pdf
  if a filename is provided.
}

\author{
  
Lukas Burger lukas.burger@fmi.ch

}

\examples{

library(MethylSeekR)

# get chromosome lengths
library("BSgenome.Hsapiens.UCSC.hg18")
sLengths=seqlengths(Hsapiens)

# read methylation data
methFname <- system.file("extdata", "Lister2009_imr90_hg18_chr22.tab",
package="MethylSeekR")
meth.gr <- readMethylome(FileName=methFname, seqLengths=sLengths)

#load CpG islands
library(rtracklayer)
session <- browserSession()
genome(session) <- "hg18"
query <- ucscTableQuery(session, table = "cpgIslandExt")
CpGislands.gr <- track(query)
genome(CpGislands.gr) <- NA
CpGislands.gr <- resize(CpGislands.gr, 5000, fix="center")

#calculate FDRs, assuming no PMDs
stats <- calculateFDRs(m=meth.gr, CGIs=CpGislands.gr)

}
