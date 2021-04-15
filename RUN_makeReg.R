rm(list=ls())
  
# Loading libraries
library(data.table)
library(GenomicRanges)
library(dplyr)
library(stringi)
library(Biostrings)

# Read source code and define output directory
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/makeRegScaffolds.R"))
source(paste0(Sys.getenv("HOME"),"/DMR-Analysis/CfromFASTAv4.R"))

#-----------------------------------------------------------------------------------
# For Arabidopsis, the FASTA files and Region files are already provided
#-----------------------------------------------------------------------------------
out.dir <-"/home/rashmi/DMR-Analysis/min.C_5/Arabidopsis/fp0.01/"
fasta.folder <-"/home/rashmi/DMR-Analysis/FASTA/"
out.name <- "Arabidopsis"
contexts <- c("CG", "CHG", "CHH","C")
makeNull <- c(TRUE, TRUE, TRUE, TRUE)
min.C <- 5
fp.rate <- 0.01


chrfiles <- list.files(fasta.folder, pattern=paste0("*.fa.gz$"), full.names = TRUE)

# I am creating a new folder 'min.C_5' here
if (!dir.exists(paste0(out.dir, "min.C_5"))) {
cat(paste0("Creating directory "))
dir.create(paste0(out.dir, "min.C_5"))
} else {
cat("directory exists!")
}

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
  ref.genome <- fread(paste0(out.dir, "min.C_5/cytosine_positions_chr", chr, ".csv", sep = ""))
  
  system.time(makeReg(ref.genome = ref.genome, 
                      contexts = contexts, 
                      makeRegnull = makeNull, 
                      chr = chr, 
                      min.C = min.C, 
                      N.boot=10^5, 
                      N.sim.C = "all", 
                      fp.rate=fp.rate, 
                      set.tol=0.01, 
                      out.dir=paste0(out.dir, "min.C_5/"), 
                      out.name=out.name))
}

#-----------------------------------------------------------------------------------
# For genomes with scaffolds
#-----------------------------------------------------------------------------------
out.dir <- "/home/rashmi/jDMR-output/min.C_5/Beech/fp0.01/"

# Either provide the longest few scaffolds. or supply all (they should be split into individual scaffolds)
# faidx --split-files scaffolds.fa

fasta.folder <- "/home/rashmi/Beech/ref_genome/fasta/"
out.name <- "Beech"
contexts <- c("CG", "CHG", "CHH")
makeNull <- c(TRUE, TRUE, TRUE)
min.C <- 5
fp.rate <- 0.01


chrfiles <- list.files(fasta.folder, pattern=paste0("*.fasta.gz$"), full.names = TRUE)

for (i in 1:length(chrfiles)) {
  tryCatch({
    fasta <-readDNAStringSet(chrfiles[i])
    #chr <- gsub("\\_size.*.fasta.gz$", "", basename(chrfiles[i]))
    chr <- gsub("\\.*.fasta.gz$", "", basename(chrfiles[i]))
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
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), chr,"\n")
    })
}
