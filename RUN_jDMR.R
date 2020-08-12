#library(devtools)
#devtools::install_github("jlab-code/jDMR")
rm(list=ls())
library(data.table)
library(dplyr)
library(methimpute)

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

#-----------------------------------------------------------------------------
# Run DMR Matrix
#-----------------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(dplyr)
library(stringr)

# makeDMRmatrix: for both population level and  pairwise control-Treatment data
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/makeDMRmatrix.R"))

input.dir <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/RegionCalls"
out.dir <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/DMRmatrix"
filelist <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/listFiles.fn"

# make binary & rc.meth.lvl matrix
makeDMRmatrix(context=c("CG"),
              chr=c("chr1","chr2","chr3","chr4","chr5"),
              samplefile=filelist,
              input.dir=input.dir,
              out.dir=out.dir)


#-----------------------------------------------------------------------------
# Run Filter DMR Matrix
#-----------------------------------------------------------------------------

rm(list=ls())
library(data.table)
library(dplyr)
library(stringr)

#--------------------------------------------
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/filterDMRmatrix.R"))
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/globFun.R"))

out.dir <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/DMRmatrix"

# for datasets with just 2 replicates (pairwise control-Treatment data) please set replicate.consensus value to 1.
filterDMRmatrix(replicate.consensus=1,
                #epiMAF.cutoff = 0.33
                epiMAF.cutoff=NULL,
                context=c("CG","CHG","CHH"),
                chr=c("chr1","chr2","chr3","chr4","chr5"),
                out.dir=out.dir)

#-----------------------------------------------------------------------------
# Run Annotate DMRs
#-----------------------------------------------------------------------------

rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(tidyr)

#Load source code
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/annotateDMRs.R", sep=""))

# annotation files
wd ="/Users/rashmi/basedir/DMRcaller"
gff.AT <- paste0(wd, "/Annotations/Arabidopsis_thaliana.TAIR10.43.gff3", sep="")
gff.TE <- paste0(wd, "/Annotations/TAIR10_TE.gff3", sep="")
gff.pr <- paste0(wd, "/Annotations/TAIR10_promoters.gff3",sep="")

#-----------------------------------------------
#myfile <- paste0(out.dir, "chr1_CG_state-calls-filtered.txt", sep="")
inputf <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/DMRmatrix"
myfiles <- list.files(inputf, pattern="*state-calls-filtered.txt", full.names = TRUE)

#-----------------------------------------------
out.dir <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/annotations"

annotateDMRs(
  gff=c(gff.AT, gff.TE, gff.pr),
  annotation=c("gene","promoters","TE"),
  file.list=myfiles,
  out.dir=out.dir)

#available annotations
#"chromosome","gene","mRNA","five_prime_UTR","exon","CDS",
#"three_prime_UTR","ncRNA_gene","lnc_RNA","miRNA","tRNA","ncRNA",
#"snoRNA","snRNA","rRNA","TE","promoters"
