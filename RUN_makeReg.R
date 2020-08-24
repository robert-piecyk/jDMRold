rm(list=ls())
  
# Loading libraries
library(data.table)
library(GenomicRanges)
library(dplyr)
library(stringi)
library(Biostrings)

# Read source code and define output directory
source("/Users/rashmi/basedir/DMRcaller/makeRegScripts/DMRs/makeReg.R")
source("/Users/rashmi/basedir/DMRcaller/makeRegScripts/DMRs/CfromFASTAv4.R")

# Edit below:
#-----------------------------------------------------------------------------------

out.dir <-"/Users/rashmi/basedir/DMRcaller/PRODUCED/min.C_5/"
fasta.folder <-"/Users/rashmi/basedir/DMRcaller/FASTA"
out.name <- "Arabidopsis"
contexts <- c("CG", "CHG", "CHH","C")
makeNull <- c(TRUE, TRUE, TRUE, TRUE)
min.C <- 5
fp.rate <- 0.01

#-----------------------------------------------------------------------------------
chrfiles <- list.files(fasta.folder, pattern=paste0("*.fa.gz$"), full.names = TRUE)
for (i in 1:length(chrfiles)) {
  fasta <-readDNAStringSet(chrfiles[i])
  chr <- gsub(".*chromosome.|\\.fa.gz$", "", basename(chrfiles[i]))
  cat(paste0("Running for chr:", chr,"\n"), sep = "")
  
  # extract cytosines from Fasta
  system.time(CfromFASTAv4(fasta = fasta, 
                           chr=chr,
                           out.dir=out.dir, 
                           write.output=TRUE))
  
  # Calling regions; calls the file created by CfromFASTAv4 
  ref.genome <- fread(paste0(out.dir, "cytosine_positions_chr", chr, ".csv", sep=""))
  system.time(makeReg(ref.genome = ref.genome, 
                      contexts = contexts, 
                      makeRegnull = makeNull, 
                      chr = chr, 
                      min.C = min.C, 
                      N.boot=10^5, 
                      N.sim.C = "all", 
                      fp.rate=fp.rate, 
                      set.tol=0.01, 
                      out.dir, 
                      out.name=out.name))
}