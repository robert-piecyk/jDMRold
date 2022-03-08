#-------------------------------------------------------------------------------
# Install jDMR
#-------------------------------------------------------------------------------
rm(list=ls())
library(devtools)
devtools::install_github("jlab-code/jDMR")
library(jDMR)

#-------------------------------------------------------------------------------
# Step 1: Run jDMR on cytosine clusters
#-------------------------------------------------------------------------------
runjDMRregions(fasta.file="/myfolder/",
               samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
               genome="Arabidopsis",
               out.dir="/myfolder/")

#-------------------------------------------------------------------------------
# Step 1: Run jDMR using grid approach
#-------------------------------------------------------------------------------
runjDMRgrid(out.dir="/myfolder/",
            fasta.file="/Annotations/TAIR10_chr_all.fa",
            samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
            min.C=10,
            genome="Arabidopsis")

#-------------------------------------------------------------------------------
# Step 2: Generate DMR matrix
#-------------------------------------------------------------------------------
makeDMRmatrix(samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
              input.dir="/myfolder/",
              out.dir="/myfolder/")

#-------------------------------------------------------------------------------
# Run this step only if you have multiple treatment groups and you want to
# perform pairwise comparisons with each of the treatment groups with control. 
# This function will split the DMR matrix into pairwise control-treatment groups. 
# Refer to section 1.2 in the manual. The metadata file requires an additional column
# "group"
#-------------------------------------------------------------------------------
split.groups(samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
             input.dir="/myfolder/",
             out.dir="/myfolder/")

#-------------------------------------------------------------------------------
# Step 3: Filter the DMR matrix
#-------------------------------------------------------------------------------
filterDMRmatrix(gridDMR=TRUE,
                data.dir="/myfolder/",
                replicate.consensus=0.8,
                samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"))

#-------------------------------------------------------------------------------
# Step 4: Generate context-specific DMRs
#-------------------------------------------------------------------------------
context.specific.DMRs(samplefiles=system.file("extdata", "listFiles2.fn", package="jDMR"),
                      data.dir="/myfolder/")

#-------------------------------------------------------------------------------
# Step 5: Annotate DMRs. Please create a new folder and move all files to be annotated 
# into the new folder
#-------------------------------------------------------------------------------
data.dir <- "/myfolder/annotate_DMRs"
gff.AT <- "/Annotations/Arabidopsis_thaliana.TAIR10.43.gff3"
gff.TE <- "/Annotations/TAIR10_TE.gff3"
gff.pr <- "/Annotations/TAIR10_promoters.gff3"

annotateDMRs(gff.files=c(gff.AT, gff.TE, gff.pr),
             annotation=c("gene","TE","promoters"), #string containing annotation types
             input.dir=data.dir,
             gff3.out=FALSE,
             out.dir=data.dir)
