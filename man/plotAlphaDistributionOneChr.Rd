\name{plotAlphaDistributionOneChr}
\alias{plotAlphaDistributionOneChr}
\title{
  Calculate and plot alpha distribution.
}

\description{

  This function calculates the alpha values for a selected chromosome
  and plots the distribution of alpha values. The shape of the
  distribution is indicative of the presence or absence of partially
  methylated domains (PMDs).

}
\usage{
plotAlphaDistributionOneChr(m, chr.sel, pdfFilename = NULL, num.cores = 1, nCGbin = 101)
}
\arguments{
  \item{m}{
    Methylation GRanges object.
  }
  \item{chr.sel}{
    Selected chromosome for which alpha values are calculated. Must be
    one of the sequence levels of m.
  }
  \item{pdfFilename}{
    Name of the pdf file in which the figure is saved. If no name is
    provided (default), the figure is printed to the screen.
  }  
  \item{num.cores}{
    The number of cores that are used for the
    calculation (default 1).
  }
  \item{nCGbin}{
    The number of CpGs in each sliding window used to calculate
    alpha (default 101). The default is highly recommended.
  }
}

\value{
  
  No return value. The function creates a figure
  showing the alpha distribution for the selected chromosome that is either
  printed to the screen (default) or saved as a pdf if a filename is
  provided.

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

#calculate alpha distribution for one chromosome
plotAlphaDistributionOneChr(m=meth.gr, chr.sel="chr22", num.cores=1)


}
