\name{segmentPMDs}
\alias{segmentPMDs}
\title{
  PMD segmenter
}
\description{
This function trains a Hidden Markov Model (HMM) to detect partially methylated domains (PMDs) in Bis-seq data.
}
\usage{
segmentPMDs(m, chr.sel, pdfFilename = NULL, seqLengths, num.cores = 1, nCGbin = 101)
}
\arguments{
  \item{m}{
    GRanges object containing the methylation data.
  }
  \item{chr.sel}{
    Chromosome on which HMM should be trained. Must be one of the
    sequence levels of m. 
  }
  \item{pdfFilename}{
    Name of the pdf file in which the figure is saved. If no name is
    provided (default), the figure is printed to the screen.
  }
  \item{seqLengths}{
    A named vector indicating the chromosome lengths of
    the genome used.
  }
  \item{num.cores}{
    The number of cores used for the calculations (default 1).
  }
  \item{nCGbin}{
    The number of CpGs in each sliding window used to calculate
    alpha (default 101). The default is highly recommended.
  }
}


\value{

  A GRanges object containing segments that partition the genome into
  PMDs and regions outside of PMDs. The object contains two metadata
  columns indicating the type of region (PMD/notPMD) and the number of
  covered (by at least 5 reads) CpGs (nCG) in the region. The function
  also creates a figure showing the inferred emission distributions of
  the HMM that is either printed to the screen (default) or saved as a
  pdf if a filename is provided.

} \author{ Lukas Burger lukas.burger@fmi.ch }

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
PMDsegments.gr <- segmentPMDs(m=meth.gr, chr.sel="chr22",
seqLengths=sLengths)


}
