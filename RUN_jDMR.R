#library(devtools)
#devtools::install_github("jlab-code/jDMR")
rm(list=ls())
library(data.table)
library(dplyr)
library(methimpute)
library(stringr)

#-----------------------------------------------------------------------------
# Step1: Run Methimpute for cytosine regions
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
filelist <- "/home/rashmi/DMR-Analysis/listFiles.fn"

#Region DMRs
##Folder containing Cytosine regions as Rdata files 
Regionsfolder <- "/home/rashmi/DMR-Analysis/min.C_5/Arabidopsis/fp0.01"
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
                  step=100,
                  genome="Arabidopsis",
                  samplefiles=filelist,
                  mincov=0,
                  nCytosines=5,
                  context=c("CG","CHG","CHH"))

#-----------------------------------------------------------------------------
# Step2: Run DMR Matrix
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
# Step3: Run Filter DMR Matrix
#-----------------------------------------------------------------------------
rm(list=ls())

source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/filterDMRmatrix.R"))
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/globFun.R"))

data.dir <- "/home/rashmi/jDMR-output"

## Please run filterDMRmatrix function based on the type of data you have.

filterDMRmatrix(replicate.consensus=NULL, # set value if your data is pairwise control-Treatment data or else set to NULL.
                # For datasets with just 2 replicates,please set replicate.consensus value to 1.
                gridDMR=TRUE, # if Region DMRs set it to FALSE else if grid DMRs set to TRUE
                epiMAF.cutoff=NULL, # set to NULL or set epiMAF.cutoff (for population data)
                data.dir=data.dir)

#-----------------------------------------------------------------------------
# Step4: Run Annotate DMRs
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
# Step 6: Methimpute to BedGraph format. 
# Converts Region level Methimpute output files to Bedgraph 
# The bedgraph outputs can be further converted to bigwig using MethylStar
#-----------------------------------------------------------------------------
rm(list=ls())
library(data.table)
 
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/MethimputeRegTobedGraph.R", sep=""))
 
wd <- "~/jDMR-output/RegionCalls"
myfiles <- list.files(wd, pattern=".txt$", full.names = TRUE)
out.dir <- "~/jDMR-output/bedgraph-out"
 
for (i in 1:length(myfiles)) {
  MethimputeRegTobedGraph.stateCalls(regfile=myfiles[i], out.dir=out.dir)
}

for (j in 1:length(myfiles)) {
  MethimputeRegTobedGraph.rcmethlvl(regfile=myfiles[j], out.dir=out.dir)
}
