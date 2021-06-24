#-------------------------------------------------------------------------------
# Install jDMR
#-------------------------------------------------------------------------------
rm(list=ls())
library(devtools)
devtools::install_github("jlab-code/jDMR")
library(jDMR)

#-------------------------------------------------------------------------------
# Create Region files for the first time (Not required if running for A.thaliana)
#-------------------------------------------------------------------------------
library(Biostrings)
library(data.table)

myoutput <- "~/jDMR-test/"

fasta <- system.file("extdata","Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa.gz", package="jDMR")
myfasta <- readDNAStringSet(fasta)
CfromFASTAv4(fasta=myfasta,
             chr=1,
             out.dir=myoutput,
             write.output=TRUE
)

ref.genome <- fread(paste0(myoutput, "/cytosine_positions_chr", 1, ".csv", sep=""))

makeReg(ref.genome=ref.genome,
        contexts=c("CG","CHG","CHH"), # all contexts can be specified as ("CG","CHG","CHH","C")
        makeRegnull=c(TRUE), # can be set to FALSE if null distribution is not required to be generated
        chr=1,
        min.C=5,
        N.boot=10^5,
        N.sim.C="all",
        fp.rate=0.01,
        set.tol=0.01,
        out.dir=myoutput,
        out.name="Arabidopsis"
)

#-------------------------------------------------------------------------------
# Call Region DMRs
#-------------------------------------------------------------------------------

myoutput <- "~/jDMR-test/"
samplefile1 <- "~/jDMR-test/listFiles.fn"
Regionsfolder <- system.file("extdata","min.C_5/fp0.01", package="jDMR")
runMethimputeRegions(Regionfiles=Regionsfolder,
                     samplefiles=samplefile1,
                     genome="Arabidopsis",
                     context=c("CG","CHG","CHH"),
                     out.dir=myoutput)

#-------------------------------------------------------------------------------
# OR, Call grid DMRs
#-------------------------------------------------------------------------------

myoutput <- "~/jDMR-test/"
samplefile1 <- "~/jDMR-test/listFiles.fn"
fasta.files <- system.file("extdata", package="jDMR")
runMethimputeGrid(fasta=fasta.files,
                  samplefiles=samplefile1,
                  genome="Arabidopsis",
                  context=c("CG","CHG","CHH"),
                  out.dir=myoutput,
                  win=100,
                  step=100,
                  mincov=0,
                  nCytosines=5)

#-------------------------------------------------------------------------------
# Make DMR matrix
#-------------------------------------------------------------------------------
makeDMRmatrix(context=c("CG","CHG","CHH"),
              samplefiles=samplefile1,
              input.dir=myoutput,
              out.dir=myoutput)

#-------------------------------------------------------------------------------
# Filter the Matrix and merge bins in case of grid DMRs
#-------------------------------------------------------------------------------

filterDMRmatrix(gridDMR=TRUE, #if region DMRs set to FALSE
                data.dir=myoutput)

##optional (by default both options set to NULL)
#replicate.consensus=8
#epiMAF.cutoff=0.33

#-------------------------------------------------------------------------------
# Annotate DMRs
#-------------------------------------------------------------------------------

# annotation files
gff.AT <- "/Annotations/Arabidopsis_thaliana.TAIR10.47.gff3"
gff.TE <- "/Annotations/TAIR10_TE.gff3"
gff.pr <- "/Annotations/TAIR10_promoters.gff3"

#Please supply the text files to be annotated in a separate folder.
#For e.g I make a new folder "mysamples". In the case of gridDMR supply the (*merged.txt) files by moving them to "mysamples" folder
mydir <- paste0(myoutput, "mysamples")

#you can specify the following available annotations. if you have your custom file let me know.
#"chromosome","gene","mRNA","five_prime_UTR","exon","CDS",
#"three_prime_UTR","ncRNA_gene","lnc_RNA","miRNA","tRNA","ncRNA",
#"snoRNA","snRNA","rRNA","TE","promoters"

annotateDMRs(gff.files=c(gff.AT, gff.TE, gff.pr),
             annotation=c("gene","promoters","TE"),
             input.dir=mydir,
             gff3.out=TRUE,
             out.dir=mydir)
