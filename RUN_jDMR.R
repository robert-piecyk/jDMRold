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
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts-local/DMRs/globFun.R"))
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts-local/DMRs/MethimputeReg.R"))
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts-local/DMRs/runMethimputeRegions.R"))

#Input files required:
#-----------------------------------------------------------------------------
##Folder containing Cytosine regions as Rdata files 
Regionsfolder <- "/Users/rashmi/basedir/DMRcaller/PRODUCED/min.C_5/fp_0.1"

##Output directory
myoutput <- "/Users/rashmi/basedir/DMRcaller/jDMR-output"

##list containing filenames with full PATH.Note that Methimpute files should have prefix "methylome" and suffix "All.txt"
filelist <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/listFiles.fn"

#-----------------------------------------------------------------------------
#Run (U,M,I) 3-state call with include.intermediate=TRUE
#probability = one of c("independent","constrained").ignore for now.

runMethimputeRegions(Regionfiles=Regionsfolder,
                     samplefiles=filelist,
                     context=c("CG","CHG","CHH"),
                     run.chr.separate=TRUE,
                     out.dir=myoutput)

#-----------------------------------------------------------------------------
# Run DMR Matrix
#-----------------------------------------------------------------------------
rm(list=ls())

# makeDMRmatrix: for both population level and  pairwise control-Treatment data
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts-local/DMRs/makeDMRmatrix.R"))

mydir <- "/Users/rashmi/basedir/DMRcaller/jDMR-output"
filelist <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/listFiles.fn"

# make binary & rc.meth.lvl matrix
makeDMRmatrix(context=c("CG","CHG","CHH"),
              samplefiles=filelist,
              input.dir=mydir,
              out.dir=mydir)

#-----------------------------------------------------------------------------
# Run Filter DMR Matrix
#-----------------------------------------------------------------------------
rm(list=ls())

source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts-local/DMRs/filterDMRmatrix.R"))
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts-local/DMRs/globFun.R"))

data.dir <- "/Users/rashmi/basedir/DMRcaller/jDMR-output"

# for datasets with just 2 replicates (pairwise control-Treatment data) please set replicate.consensus value to 1.
# either specify value or set replicate.consensus=NULL or epiMAF.cutoff=NULL
filterDMRmatrix(replicate.consensus=1,
                #epiMAF.cutoff = 0.33,
                epiMAF.cutoff=NULL,
                data.dir=data.dir)

#-----------------------------------------------------------------------------
# Run Annotate DMRs
#-----------------------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(rtracklayer)
library(tidyr)

#Load source code
source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/annotateDMRs.R", sep=""))

# gff3 annotation files. Supply as one gff3 file
wd ="/Users/rashmi/basedir/DMRcaller"
gff.AT <- paste0(wd, "/Annotations/Arabidopsis_thaliana.TAIR10.43.gff3", sep="")
gff.TE <- paste0(wd, "/Annotations/TAIR10_TE.gff3", sep="")
gff.pr <- paste0(wd, "/Annotations/TAIR10_promoters.gff3",sep="")
gff <- c(gff.AT, gff.TE, gff.pr)
input.gff <- lapply(gff, function(x){ 
  import.gff3(x, colnames=c("type", "ID")) 
})
merged.gff <- do.call(c, input.gff)

# Check Annotation levels here. Supply annotation terms for e.g genes, TEs
levels(elementMetadata(merged.gff)[,"type"])
#available annotations
#"chromosome","gene","mRNA","five_prime_UTR","exon","CDS",
#"three_prime_UTR","ncRNA_gene","lnc_RNA","miRNA","tRNA","ncRNA",
#"snoRNA","snRNA","rRNA","TE","promoters"

# Path to filtered DMR files with 3 columns: seqnames, start, end
inputf <- "/Users/rashmi/basedir/DMRcaller/CytosineRegions_background/FP0.1"
myfiles <- list.files(inputf, pattern="*.txt", full.names = TRUE)

# Output Path
out.dir <- "/Users/rashmi/basedir/DMRcaller/CytosineRegions_background/FP0.1"

annotateDMRs(gff=merged.gff,
             annotation=c("gene","promoters","TE"),
             file.list=myfiles,
             gff3.out=TRUE,
             out.dir=out.dir)

#-----------------------------------------------------------------------------
# Methimpute to BedGraph format. This is only for Region level output files
# The bedgraph outputs can be further converted to bigwig using MethylStar
#-----------------------------------------------------------------------------
rm(list=ls())
library(data.table)

source(paste0(Sys.getenv("HOME"),"/basedir/DMRcaller/makeRegScripts/DMRs/MethimputeRegTobedGraph.R", sep=""))

wd <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/RegionCalls"
myfiles <- list.files(wd, pattern=".txt$", full.names = TRUE)
out.dir <- "/Users/rashmi/basedir/DMRcaller/jDMR-output/bedgraph"

for (i in 1:length(myfiles)) {
  MethimputeRegTobedGraph(regfile=myfiles[i],
                          out.dir=out.dir)
}

