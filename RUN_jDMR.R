#library(devtools)
#devtools::install_github("jlab-code/jDMR")
rm(list=ls())
library(data.table)
library(dplyr)
library(methimpute)
library(stringr)

#-----------------------------------------------------------------------------
# Run Methimpute for cytosine regions
#-----------------------------------------------------------------------------
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/globFun.R"))
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/MethimputeReg.R"))
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/runMethimpute.R"))

#Input files: 
#-----------------------------------------------------------------------------
##Output directory
myoutput <- "/home/rashmi/jDMR-output"

##list containing filenames with full PATH.
##Note that Methimpute files should have prefix "methylome" and suffix "All.txt"
filelist <- "/home/rashmi/DMR-Analysis/listFiles-AT.fn"

#-----------------------------------------------------------------------------
#Run (U,M,I) 3-state call with include.intermediate=TRUE
#probability = one of c("independent","constrained").ignore for now.

#Region DMRs
##Folder containing Cytosine regions as Rdata files 
Regionsfolder <- "/home/rashmi/jDMR-Analysis/min.C_5/fp_0.1"
runMethimputeRegions(Regionfiles=Regionsfolder,
                     samplefiles=filelist,
                     genome="Arabidopsis",
                     context=c("CG","CHG","CHH"),
                     out.dir=myoutput)

#grid DMRs
#Here you need to specify the chromosome lengths. For Arabidopsis you can use the ones provided.
chrs <- c(chr1=30427671, chr2=19698289, chr3=23459830, chr4=18585056, chr5=26975502)
runMethimputeGrid(out.dir=myoutput, 
                  chrfile=chrs,
                  win=100, 
                  step=50,
                  genome="Arabidopsis",
                  samplefiles=filelist,
                  mincov=10,
                  nCytosines=10,
                  context=c("CG","CHG","CHH"))

#-----------------------------------------------------------------------------
# Run DMR Matrix
#-----------------------------------------------------------------------------
rm(list=ls())

# makeDMRmatrix: for both population level and  pairwise control-Treatment data
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/makeDMRmatrix.R"))

mydir <- "/home/rashmi/jDMR-output"
filelist <- "/home/rashmi/DMR-Analysis/listFiles-AT.fn"

# make binary & rc.meth.lvl matrix
makeDMRmatrix(context=c("CG","CHG","CHH"),
              samplefiles=filelist,
              input.dir=mydir,
              out.dir=mydir)

#-----------------------------------------------------------------------------
# Run Filter DMR Matrix
#-----------------------------------------------------------------------------
rm(list=ls())

source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/filterDMRmatrix.R"))
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/globFun.R"))

data.dir <- "/home/rashmi/jDMR-output"

# either specify value or set replicate.consensus=NULL or epiMAF.cutoff=NULL. Please run filterDMRmatrix function based on the type of data you have.

# 1) Region DMRs for Population data
filterDMRmatrix(replicate.consensus=NULL,
                gridDMR=FALSE,
                #epiMAF.cutoff is only for population data
                epiMAF.cutoff = 0.33,
                data.dir=data.dir)

# 2) Region DMRs for pairwise control-Treatment data
# for datasets with just 2 replicates (pairwise control-Treatment data) please set replicate.consensus value to 1.
filterDMRmatrix(replicate.consensus=1,
                gridDMR=FALSE,
                epiMAF.cutoff=NULL,
                data.dir=data.dir)

# 3) grid DMRs for Population data
filterDMRmatrix(replicate.consensus=NULL,
                gridDMR=TRUE,
                #epiMAF.cutoff is only for population data
                epiMAF.cutoff=0.33,
                data.dir=data.dir)

# 4) grid DMRs for pairwise control-Treatment data
filterDMRmatrix(replicate.consensus=1,
                gridDMR=TRUE,
                epiMAF.cutoff=NULL,
                data.dir=data.dir)


#-----------------------------------------------------------------------------
# Run Annotate DMRs
#-----------------------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)
library(dplyr)

#Load source code
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/annotateDMRs.R", sep=""))

wd ="/home/rashmi"
#Please supply the text files to be annotated in a separate folder. For e.g I make a new folder "mysamples". In the case of gridDMR supply the (*merged.txt) files by moving them to "mysamples" folder
input.dir <- "/home/rashmi/jDMR-output/mysamples"
out.dir <- "/home/rashmi/jDMR-output/mysamples"

gff.AT <- paste0(wd, "/Annotations/Arabidopsis_thaliana.TAIR10.47.gff3", sep="")
gff.TE <- paste0(wd, "/Annotations/TAIR10_TE.gff3", sep="")
gff.pr <- paste0(wd, "/Annotations/TAIR10_promoters.gff3",sep="")

#you can specify the following available annotations. if you have your custom file let me know.
#"chromosome","gene","mRNA","five_prime_UTR","exon","CDS",
#"three_prime_UTR","ncRNA_gene","lnc_RNA","miRNA","tRNA","ncRNA",
#"snoRNA","snRNA","rRNA","TE","promoters"

annotateDMRs(gff.files=c(gff.AT, gff.TE, gff.pr),
             annotation=c("gene","promoters","TE"),
             input.dir=input.dir,
             gff3.out=TRUE,
             out.dir=out.dir)

#-----------------------------------------------------------------------------
# Methimpute to BedGraph format. This is only for Region level output files
# The bedgraph outputs can be further converted to bigwig using MethylStar
#-----------------------------------------------------------------------------
# rm(list=ls())
# library(data.table)
# 
# source(paste0(Sys.getenv("HOME"),"/basedir/jDMR-scripts/MethimputeRegTobedGraph.R", sep=""))
# 
# wd <- "~/basedir/jDMR-output/100win_50Stepsize/all-samples/RegionCalls"
# myfiles <- list.files(wd, pattern=".txt$", full.names = TRUE)
# out.dir <- "~/basedir/jDMR-output/100win_50Stepsize/all-samples/bedgraph"
# 
# for (i in 1:length(myfiles)) {
#   MethimputeRegTobedGraph.stateCalls(regfile=myfiles[i],
#                           out.dir=out.dir)
# }

