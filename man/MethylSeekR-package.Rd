\name{MethylSeekR-package}
\alias{MethylSeekR-package}
\alias{MethylSeekR}
\docType{package}
\title{
Segmentation of Bis-seq methylation data
}
\description{
This is a package for the discovery of regulatory regions from Bis-seq data
}
\details{
\tabular{ll}{
Package: \tab MethylSeekR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-10-10\cr
License: \tab None\cr
}


}
\author{
Lukas Burger

Maintainer: Lukas Burger <Lukas.Burger@fmi.ch>
}
\references{
Stadler, Murr, Burger et al, DNA-binding factors shape the mouse
methylome at distal regulatory regions, Nature 2011. 
}
\keyword{ package }
\examples{

library(MethylSeekR)

# get chromosome lengths
library("BSgenome.Hsapiens.UCSC.hg38")
sLengths=seqlengths(Hsapiens)

# read methylation data
methFname <- system.file("extdata", "Lister2009_imr90_hg38_chr22.tab",
package="MethylSeekR")
meth.gr <- readMethylome(FileName=methFname, seqLengths=sLengths)

#read SNP data
snpFname <- system.file("extdata", "SNVs_hg38_chr22.tab",
package="MethylSeekR")
snps.gr <- readSNPTable(FileName=snpFname, seqLengths=sLengths)

# remove SNPs
meth.gr <- removeSNPs(meth.gr, snps.gr)

#calculate alpha distribution for one chromosome
plotAlphaDistributionOneChr(m=meth.gr, chr.sel="chr22", num.cores=1)

#segment PMDs
PMDsegments.gr <- segmentPMDs(m=meth.gr, chr.sel="chr22", num.cores=1,
seqLengths=sLengths)

#plot PMD segmentation examples
plotPMDSegmentation(m=meth.gr, segs=PMDsegments.gr, numRegions=1)

#save PMD segments
savePMDSegments(PMDs=PMDsegments.gr, GRangesFilename="PMDs.gr.rds",
 TableFilename="PMDs.tab")

#load CpG islands
library(rtracklayer)
session <- browserSession()
genome(session) <- "hg38"
query <- ucscTableQuery(session, table = "cpgIslandExt")
CpGislands.gr <- track(query)
genome(CpGislands.gr) <- NA
CpGislands.gr <- resize(CpGislands.gr, 5000, fix="center")

#calculate FDRs
stats <- calculateFDRs(m=meth.gr, CGIs=CpGislands.gr,
PMDs=PMDsegments.gr, num.cores=1)

# select FDR cut-off and determine segmentation parameters
FDR.cutoff <- 5 
m.sel <- 0.5 
n.sel=as.integer(names(stats$FDRs[as.character(m.sel), ]
[stats$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])

#segment UMRs and LMRs 
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel,
nCpG.cutoff=n.sel, PMDs=PMDsegments.gr, num.cores=1,
myGenomeSeq=Hsapiens, seqLengths=sLengths)

#plot final segmentation including PMDs, UMRs and LMRs
plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr, PMDs=PMDsegments.gr,
meth.cutoff=m.sel, numRegions=1)

#save UMRs and LMRs
saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename="UMRsLMRs.gr.rds",
 TableFilename="UMRsLMRs.tab")


}
