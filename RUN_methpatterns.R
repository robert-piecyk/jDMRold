rm(list=ls())
library(data.table)
library(plyr)
library(dplyr)
library(mgcv)

source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/methpatterns.R"))

#-----------------------------------------------------------------------------
# Extract patterns
#-----------------------------------------------------------------------------

##Output directory containing methimpute region level calls
methout <- "/Users/rashmi/basedir/DMRcaller/test/chr1_CG"
methpatterns(methout=methout, 
             context=c("CG","CHG","CHH"),
             #chr <- c("chr1"),
             chr <- c("chr1","chr2","chr3","chr4","chr5"),
             out.dir=methout, 
             WT="SRR534177")
