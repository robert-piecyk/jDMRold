#' Extract cytosines from FASTA
#'
#' This function slides along the genome and extracts all Cytosines
#' @param fasta
#' @param chr
#' @param out.dir
#' @param write.output
#' @import Biostrings
#' @import utils
#' @import data.table
#' @export
#'

### CfromFASTAv1:
### - Option to return tri-nucleotide contexts if needed
### - Bioconstrings for finding positions of contexts

 CfromFASTAv4<-function(fasta, chr, out.dir, write.output)
 {

   # Defining function for cleaning objects in the global environment from within the function
   CleanEnvir <- function(pattern)
   {
     rm(list = ls(envir=globalenv())[
       grep(pattern, ls(envir=globalenv()))], envir = globalenv())
   }

      cat("- Processing chr:", chr, "\n")
      cat("-------------------------", "\n")

      cat("- Converting DNA .....", "\n")
      fasta.minus<-Biostrings::complement(fasta)
      fasta.plus<-Biostrings::DNAString(paste(fasta, collapse=""))
      fasta.minus<-Biostrings::DNAString(paste(fasta.minus, collapse=""))

      # Removing the .fasta file from global environment
      CleanEnvir(pattern="fasta")

      # Listing the cytosine patterns by context
      plus.pattern<-list("CG",  c("CAG", "CTG", "CCG"), c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC"))
      names(plus.pattern)<-c("CG", "CHG", "CHH")

      minus.pattern<-list("GC", c("GTC","GAC", "GCC"), c("AAC", "TAC", "CAC", "ATC", "TTC", "CTC", "ACC", "TCC","CCC"))
      names(minus.pattern)<-c("CG", "CHG", "CHH")

      cytosine.patterns<-list(plus.pattern, minus.pattern)
      names(cytosine.patterns)<-c("+", "-")

        collect.positions<-list()
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
                     match.out<-data.table(rep(chr, length(match.out)), match.out, rep(names(cytosine.patterns)[[s]], length(match.out)),
                                           rep(names(strand.temp)[c], length(match.out)))
                     colnames(match.out)<-c("chr", "pos", "strand", "context")
                     collect.positions[[counter]]<-match.out
                   }
                   if (names(cytosine.patterns)[[s]] == "-")
                   {
                     match.out<-Biostrings::matchPattern(context.temp[p], subject=fasta.minus)
                     match.out<-as.data.frame(slot(match.out, name="ranges"))[,2]
                     match.out<-data.table(rep(chr, length(match.out)), match.out, rep(names(cytosine.patterns)[[s]], length(match.out)),
                                           rep(names(strand.temp)[c], length(match.out)))
                     colnames(match.out)<-c("chr", "pos", "strand", "context")
                     collect.positions[[counter]]<-match.out
                   }
                }

              }

        }

        cat("- Combining contexts .....", "\n")
        C.all<-data.table(do.call("rbind", collect.positions))

        # Reading out the data
        if (write.output == TRUE)
        {
          cat("- Writing out file .....", "\n")
          fwrite(C.all, paste(out.dir, "/cytosine_positions_chr", chr, ".csv", sep=""), row.names = FALSE)
        }

        #return(C.all)

        # Cleaning up
        rm(fasta, fasta.plus, fasta.minus, collect.positions, match.out, C.all)


 }




