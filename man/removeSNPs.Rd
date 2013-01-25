\name{removeSNPs}
\alias{removeSNPs}
\title{
  Remove CpGs overlapping SNPs
}
\description{
 Removes CpGs that overlap with SNPs from methylation GRanges object.
}
\usage{
removeSNPs(m, snps)
}
\arguments{
  \item{m}{
    GRanges object containing the methylation data.
}
  \item{snps}{
    GRanges object containing the SNPs.
}
}
\value{
The methylation GRanges object (m) with all CpGs overlapping SNPs removed.
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

#read SNP data
snpFname <- system.file("extdata", "SNVs_hg18_chr22.tab",
package="MethylSeekR")
snps.gr <- readSNPTable(FileName=snpFname, seqLengths=sLengths)

# remove SNPs
meth.gr <- removeSNPs(meth.gr, snps.gr)


}


