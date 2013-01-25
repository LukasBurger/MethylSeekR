\name{readMethylome}
\alias{readMethylome}
\title{
Load Bis-seq data
}
\description{
Loading Bis-seq data from tab-delimited file or saved GRanges object
}
\usage{
readMethylome(FileName, seqLengths, format = "text")
}

\arguments{


  \item{FileName}{File name.}

  \item{seqLengths}{A named vector indicating the chromosome lengths of
    the genome used.}

  \item{format}{File format. If format is set to "text" (default), the
  argument FileName should refer to a tab-delimited text file in the
  format: chromosome position T M, where each line stands for a CpG, the
  position refers to the C of the CpG (on the plus strand), T is the
  total number of reads (total counts) covering the CpG and M the total
  number of reads without C to T conversion at the C of the CpG
  (methylation counts). If format="GRanges", the file is assumed to be a
  GRanges object, containing T and M as first and second data-value
  entries, saved in rds format.}

}

\value{
  A GRanges object containing the coordinates, total (T) and methylated
  counts (M)
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


}
