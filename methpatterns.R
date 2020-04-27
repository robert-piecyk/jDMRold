
myfunc <- function(x) {
  as.integer(sum(x * 2^(rev(seq_along(x)) - 1)))
}

diffr <- function(x, i=1){
  n = length(x)
  dff = x[i+1:n] - x[i]
  return(dff)
}

methpatterns <- function(methout, context, chr, out.dir, WT) {
  for (i in  1:length(context)){
    for (j in 1:length(chr)){
      pattern=paste0("_1_All_", chr[j], "_", context[i] ,".txt")
      cat(paste0("Running for ", chr[j], "_", context[i], "\n"), sep = "")
      myfiles <- list.files(methout, pattern=pattern, full.names = TRUE)
      
      filelist <- lapply(myfiles, function(i){
        f <- fread(i, header=FALSE, skip=1)
        fname <- gsub(pattern = pattern, "", basename(i))
        f.df <- f[,c("V1","V2","V3","V7")]
        colnames(f.df)[4] <- fname
        return(f.df)
      })
      
      df1 <- Reduce(function(x, y) {
        dplyr::inner_join(x, y, by=c("V1","V2","V3"))
      } , filelist)
      
      
      #arranging the columns
      df_WT <- df1 %>% dplyr::select("V1","V2","V3", contains(WT), everything())
      
      #Replace data.frame with -1, 1 depending on the condition
      mymat <- df_WT
      for (k in 1:nrow(mymat)) {
        for (l in 5:ncol(mymat)){
          if (mymat[k,l] < mymat[k,4]){
            mymat[k,l]=-1
          } else if (mymat[k,l] > mymat[k,4]){
            mymat[k,l]=1
          }
        }
      }
      m <- ifelse(mymat[,5:ncol(mymat)] == -1, yes = 0, no = 1)
      mm <- uniquecombs(m)
      
      #mx <- apply(mm, MARGIN = 1, myfunc)
      mx <- apply(mm[,c(1:ncol(mm))], 1, paste, collapse="")
      mxVec <- mx[attr(mm, "index")] 
      
      #get patterns
      pat <- cbind(mm, mx)
      pat <- data.frame(pat)
      #pat$pattern <- apply(pat[,c(1:ncol(pat)-1)], 1, paste, collapse="")
      #pat <- pat %>% select("mx","pattern")
      
      ndf_WT <- cbind(df_WT, mxVec)
      colnames(ndf_WT)[ncol(ndf_WT)] <- "Pattern"
      
      tot <- nrow(mymat)
      patternCounts <- as.data.frame(table(as.factor(mxVec)), stringsAsFactors = FALSE)
      patternCounts$density <- patternCounts$Freq/tot
      colnames(patternCounts)[1] <- "Pattern"
      colnames(patternCounts)[2] <- paste0(context[i],"-","Freq")
      colnames(patternCounts)[3] <- paste0(context[i],"-","density")
      #final.df <- merge(patternCounts, pat, by.x="Pattern.int", by.y="mx")
      #writing out the matrix
      fwrite(x=mymat, file=paste0(out.dir, "/", chr[j], "_", context[i],"_mat.txt"), quote=FALSE, 
             row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=ndf_WT, file=paste0(out.dir, "/", chr[j], "_", context[i],"_vals.txt"), quote=FALSE, 
             row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=patternCounts, file=paste0(out.dir, "/", chr[j], "_", context[i],"_meth-patterns-freq.txt"), 
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    }
  }
}