rm(list=ls())

library(data.table)
library(dplyr)
library(rtracklayer)
library(stringr)

# makeDMRmatrix: for both population level and  pairwise control-Treatment data
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/makeDMRmatrix.R"))

context <- c("CG","CHG","CHH")
chr <- c("chr1","chr2","chr3","chr4","chr5")

# Filelist can be easily created in linux by typing : ls -d $PWD/whateverfolder/*.txt
samplefile <- "/Users/rashmi/basedir/DMRcaller/test/listFiles.fn"
out.dir <- "/Users/rashmi/basedir/DMRcaller/test/DMRs/"

# make binary & rc.meth.lvl matrix
makeDMRmatrix(context=context,
              chr=chr,
              samplefile=samplefile,
              out.dir=out.dir)

#--------------------------------------------
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/filterDMRmatrix.R"))

# filterDMRmatrix : only for pairwise control-Treatment data
context <- c("CG","CHG","CHH")
chr <- c("chr1","chr2","chr3","chr4","chr5")
out.dir <- "/Users/rashmi/basedir/DMRcaller/test/DMRs/"
filterDMRmatrix(context=context,
         chr=chr,
         out.dir=out.dir)

