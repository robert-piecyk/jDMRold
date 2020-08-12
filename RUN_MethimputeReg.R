#library(devtools)
#devtools::install_github("jlab-code/jDMR")
rm(list=ls())
library(data.table)
library(plyr)
library(dplyr)
library(methimpute)
library(stringr)

runMethimputeRegions <- function(Methfiles, Regionfiles, context, fit.plot, out.dir) {
  filelist <- fread(Methfiles, header=TRUE)
  for (i in 1:length(filelist$file)){
    for (j in 1:length(context)){
      methfn <- gsub(".*methylome_|\\_All.txt$", "", filelist$file[i])
      Regfiles <- list.files(Regionfiles, pattern=paste0("_", context[j], ".Rdata"), full.names = TRUE)
      for (k in 1:length(Regfiles)){
        
        tmp <- gsub(".*Arabidopsis_regions_|\\.Rdata$", "", Regfiles[k])
        chr <- gsub(paste0("_", context[j]), "", tmp)
        
        name <- paste0(methfn, "_", chr)
        cat(paste0("Running file ", methfn, " for context ", context[j], " and ", chr,"\n"), sep = "")
        makeMethimpute(
          df=filelist$file[i],
          context=context[j],
          refRegion=Regfiles[k],
          fit.plot=fit.plot,
          include.intermediate=FALSE, 
          probability="constrained",
          out.dir=out.dir,
          fit.name=paste0(methfn, "_", context[j], "_", chr),
          name=name)
      }
    }
  }
}

#-----------------------------------------------------------------------------
# Run Methimpute for cytosine regions
#-----------------------------------------------------------------------------
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/globFun.R"))
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/MethimputeReg.R"))

#Input files required:
#-----------------------------------------------------------------------------
##Folder containing Cytosine regions as Rdata files 
Regionsfolder <- "/Users/rashmi/basedir/DMRcaller/PRODUCED/min.C_5/fp_0.1"

##Output directory
myoutput <- "/Users/rashmi/basedir/DMRcaller/jDMR-output"

##list containing filenames with full PATH
filelist <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/listFiles.fn"

#-----------------------------------------------------------------------------
#Run (U,M,I) 3-state call with include.intermediate=TRUE
#probability = one of c("independent","constrained").ignore for now.

runMethimputeRegions(Regionfiles=Regionsfolder,
                     Methfiles=filelist,
                     context=c("CG"),
                     fit.plot=FALSE,
                     out.dir=myoutput)

