\name{savePMDSegments}
\alias{savePMDSegments}

\title{Save PMD segments}

\description{Save PMD segments in rds format and as tab-delimited file.}

\usage{savePMDSegments(PMDs, GRangesFilename, TableFilename)}

\arguments{

  \item{PMDs}{GRanges object containing the PMD segmentation. Return
    value of the segmentPMDs function (see example).}

  \item{GRangesFilename}{Filename of the GRanges object.}

  \item{TableFilename}{Filename of the PMD table.}

}

\value{
  No return value.
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

#segment PMDs
PMDsegments.gr <- segmentPMDs(m=meth.gr, chr.sel="chr22",
seqLengths=sLengths)

#save PMD segments
savePMDSegments(PMDs=PMDsegments.gr, GRangesFilename="PMDs.gr.rds",
 TableFilename="PMDs.tab")


}
