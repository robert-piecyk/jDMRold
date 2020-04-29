#library(devtools)
#devtools::install_github("jlab-code/jDMR")
rm(list=ls())
library(data.table)
library(plyr)
library(dplyr)
library(methimpute)
library(stringr)

#-----------------------------------------------------------------------------
# Run Methimpute for cytosine regions
#-----------------------------------------------------------------------------
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/globFun.R"))
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/MethimputeReg.R"))

#Input files required:
#-----------------------------------------------------------------------------
##Folder containing Cytosine regions as Rdata files 
Regionsfolder <- "/Users/rashmi/basedir/DMRcaller/PRODUCED/min.C_6"
##Output directory
myoutput <- "/Users/rashmi/basedir/DMRcaller/"
##Folder containing Methimpute outputs 
Methfiles <- "/Users/rashmi/basedir/DMRcaller/methimpute-out"
#-----------------------------------------------------------------------------
#Run (U,M,I) 3-state call with include.intermediate=TRUE
#probability = one of c("independent","constrained").

runMethimputeRegions(Regionfiles=Regionsfolder,
                     Methfiles=Methfiles,
                     context=c("CG","CHG","CHH"),
                     fit.plot=FALSE,
                     include.intermediate=FALSE,
                     probability="constrained",
                     out.dir=myoutput)

