#-------------------------------------------------------------------------------
# Install jDMR
#-------------------------------------------------------------------------------
rm(list=ls())
library(devtools)
devtools::install_github("jlab-code/jDMR")
library(jDMR)

#-------------------------------------------------------------------------------
# Run jDMR on cytosine clusters
#-------------------------------------------------------------------------------
runjDMRregions(fasta.file="/myfolder/",
               samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
               genome="Arabidopsis",
               out.dir="/myfolder/")

#-------------------------------------------------------------------------------
# Run jDMR using grid approach
#-------------------------------------------------------------------------------
runjDMRgrid(out.dir="/myfolder/",
            fasta.file="/Annotations/TAIR10_chr_all.fa",
            samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
            min.C=10,
            genome="Arabidopsis")

#-------------------------------------------------------------------------------
# Generate DMR matrix
#-------------------------------------------------------------------------------
makeDMRmatrix(samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
              input.dir="/myfolder/",
              out.dir="/myfolder/")

#-------------------------------------------------------------------------------
# Run this step only if you have multiple treatment groups and you want to
# perform pairwise comparisons with control. This function will split the DMR matrix
# into pairwise control-treatment groups. Don't forget to edit the metadata file before
# running the step
#-------------------------------------------------------------------------------
split.groups(samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
             input.dir="/myfolder/",
             out.dir="/myfolder/")

#-------------------------------------------------------------------------------
# Filter the DMR matrix
#-------------------------------------------------------------------------------
filterDMRmatrix(gridDMR=TRUE,
                data.dir="/myfolder/",
                replicate.consensus=0.8,
                samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"))

#-------------------------------------------------------------------------------
# Generate context-specific DMRs
#-------------------------------------------------------------------------------
context.specific.DMRs(samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
                      data.dir="/myfolder/")

#-------------------------------------------------------------------------------
# Annotate DMRs
#-------------------------------------------------------------------------------
data.dir <- "/myfolder/"
gff.AT <- "/Annotations/Arabidopsis_thaliana.TAIR10.43.gff3"
gff.TE <- "/Annotations/TAIR10_TE.gff3"
gff.pr <- "/Annotations/TAIR10_promoters.gff3"

annotateDMRs(gff.files=c(gff.AT, gff.TE, gff.pr),
             annotation=c("gene","promoters","TE"),
             input.dir=data.dir,
             gff3.out=FALSE,
             out.dir=data.dir)
