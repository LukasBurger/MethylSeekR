\name{readSNPTable}
\alias{readSNPTable}
\title{
Load SNP table
}
\description{
Loading SNPs from tab-delimited file or saved GRanges object
}
\usage{
readSNPTable(FileName, seqLengths, format = "text")
}


\arguments{

  \item{FileName}{File Name.}

  \item{seqLengths}{A named vector indicating the chromosome lengths of
  the genome used.}

  \item{format}{File format. If format is set to "text", the argument
  FileName should refer to a tab-delimited text file in the format:
  chromosome position, where each line stands for a SNP. If
  format="GRanges", the file is assumed to be a GRanges object, containing
  the SNP coordinates, saved in rds format.}

}


\value{
  A GRanges object containing the coordinates of the SNPs.
}
\author{
  Lukas Burger lukas.burger@fmi.ch
}

\examples{

library(MethylSeekR)

# get chromosome lengths
library("BSgenome.Hsapiens.UCSC.hg38")
sLengths=seqlengths(Hsapiens)

#read SNP data
snpFname <- system.file("extdata", "SNVs_hg38_chr22.tab",
package="MethylSeekR")
snps.gr <- readSNPTable(FileName=snpFname, seqLengths=sLengths)


}
