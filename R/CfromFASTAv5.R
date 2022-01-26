#rm(list=ls())
# library(Biostrings)
### CfromFASTAv1:
### - Option to return tri-nucleotide contexts if needed
### - Bioconstrings for finding positions of contexts
#fasta<-"/Users/rashmi/basedir/DMRcaller/Annotations/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"

CfromFASTAv5<-function(fasta, out.dir=NULL, write.output=NULL){

  # Defining function for cleaning objects in the global environment from within the function
  #CleanEnvir <- function(pattern)
  #{
  #  rm(list = ls(envir=globalenv())[
  #    grep(pattern, ls(envir=globalenv()))], envir = globalenv())
  #}

  all.chr <- list()
  myfasta <- Biostrings::readDNAStringSet(fasta)

  # Listing the cytosine patterns by context
  plus.pattern<-list("CG",  c("CAG", "CTG", "CCG"), c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC"))
  names(plus.pattern)<-c("CG", "CHG", "CHH")

  minus.pattern<-list("GC", c("GTC","GAC", "GCC"), c("AAC", "TAC", "CAC", "ATC", "TTC", "CTC", "ACC", "TCC","CCC"))
  names(minus.pattern)<-c("CG", "CHG", "CHH")

  cytosine.patterns<-list(plus.pattern, minus.pattern)
  names(cytosine.patterns)<-c("+", "-")

  for (chr in seq_along(myfasta)){
    mychr <- names(myfasta[chr])

    cat("- Processing chr:", mychr, "\n")
    cat("-------------------------", "\n")

    cat("- Converting DNA .....", "\n")
    #fasta.minus <- Biostrings::complement(myfasta[chr])
    fasta.plus <- Biostrings::DNAString(paste(myfasta[chr], collapse=""))
    fasta.minus <- Biostrings::DNAString(paste(Biostrings::complement(myfasta[chr]), collapse=""))

    # Removing the .fasta file from global environment
    #CleanEnvir(pattern="fasta")

    collect.positions <- list()
    counter=0

    for (s in 1:length(cytosine.patterns))
    {
      strand.temp<-cytosine.patterns[[s]]

      for (c in 1:length(strand.temp))
      {
        cat("- Scanning", names(strand.temp)[c], names(cytosine.patterns)[s], "strand", ".....", "\n")
        context.temp<-strand.temp[[c]]

        for (p in 1:length(context.temp))
        {
          counter<-counter+1

          if (names(cytosine.patterns)[[s]] == "+")
          {
            match.out<-Biostrings::matchPattern(context.temp[p], subject=fasta.plus)
            match.out<-as.data.frame(slot(match.out, name="ranges"))[,1]
            match.out<-data.table(rep(mychr, length(match.out)), match.out, rep(names(cytosine.patterns)[[s]], length(match.out)),
                                  rep(names(strand.temp)[c], length(match.out)))
            colnames(match.out) <- c("chr", "pos", "strand", "context")
            collect.positions[[counter]]<-match.out
          }
          if (names(cytosine.patterns)[[s]] == "-")
          {
            match.out<-Biostrings::matchPattern(context.temp[p], subject=fasta.minus)
            match.out<-as.data.frame(slot(match.out, name="ranges"))[,2]
            match.out<-data.table(rep(mychr, length(match.out)), match.out, rep(names(cytosine.patterns)[[s]], length(match.out)),
                                  rep(names(strand.temp)[c], length(match.out)))
            colnames(match.out) <- c("chr", "pos", "strand", "context")
            collect.positions[[counter]]<-match.out
          }
        }
      }
    }
    cat("- Combining contexts .....", "\n")
    C.all <- data.table(do.call("rbind", collect.positions))
    all.chr[[chr]] <- C.all

    # Cleaning up
    rm(fasta.plus, fasta.minus, collect.positions, match.out, C.all)
  }
  out.C <- rbindlist(all.chr)
  return(out.C)
}
