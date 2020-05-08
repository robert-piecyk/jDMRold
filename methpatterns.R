
myfunc <- function(x) {
  as.integer(sum(x * 2^(rev(seq_along(x)) - 1)))
}

diffr <- function(x, i=1){
  n = length(x)
  dff = x[i+1:n] - x[i]
  return(dff)
}

methpatterns <- function(methout, context, chr, out.dir, WT) {
  
  val.mat <- list()
  bin.mat <- list()
  freq.tbl <- list()
  tot <- 0
  
  for (i in  1:length(context)){
    cat(paste0("Running for context ", context[i], " .........\n"), sep="")
    cat("\n")
    
    for (j in 1:length(chr)){
      pattern=paste0("_1_All_", chr[j], "_", context[i] ,".txt")
      cat(paste0("Merging samples into one dataframe for ", chr[j] , "\n"), sep = "")
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
      
      df_WT.1 <- cbind(df_WT, mxVec)
      mymat.1 <- cbind(mymat, mxVec)
      colnames(df_WT.1)[ncol(df_WT.1)] <- "Pattern"
      colnames(mymat.1)[ncol(mymat.1)] <- "Pattern"
      
      tot <- tot + nrow(mymat[i])
      
      patternCounts <- as.data.frame(table(as.factor(mxVec)), stringsAsFactors = FALSE)
      
      colnames(patternCounts)[1] <- "Pattern"
      
      val.mat[[j]] <- data.frame(df_WT.1)
      bin.mat[[j]] <- data.frame(mymat.1)
      freq.tbl[[j]] <- data.frame(patternCounts)
      
      #final.df <- merge(patternCounts, pat, by.x="Pattern.int", by.y="mx")
    }
    
    cat("Now, merging dataframes for each chromosome into one file........\n")
    cat("\n")
    
    print (paste0("Total regions/rows for all chromosomes = ", tot, sep=""))
    
    val.mat.final <- rbindlist(val.mat)
    bin.mat.final <- rbindlist(bin.mat)
    
    freq.tbl.final <- rbindlist(freq.tbl)
    freq.tbl.final.1 <- aggregate(. ~ Pattern, data=freq.tbl.final, FUN=sum)
    freq.tbl.final.1$density <- freq.tbl.final.1$Freq/tot
    
    #writing out the matrix
    fwrite(x=val.mat.final, file=paste0(out.dir, "/", context[i],"_All_vals.txt"), quote=FALSE, 
           row.names=FALSE, col.names=TRUE, sep="\t")
    fwrite(x=bin.mat.final, file=paste0(out.dir, "/", context[i],"_All_mat.txt"), quote=FALSE, 
           row.names=FALSE, col.names=TRUE, sep="\t")
    fwrite(x=freq.tbl.final.1, file=paste0(out.dir, "/", context[i],"_All_methpatterns-freq.txt"), quote=FALSE, 
           row.names=FALSE, col.names=TRUE, sep="\t")
  }
}


filter.methpatterns <- function(val.matrix, freq, density.cutoff, out.dir){
  f1 <- fread(freq, header=TRUE, colClasses=c("Pattern"="character"))
  f2 <- fread(val.matrix, header=TRUE, colClasses=c("Pattern"="character"), stringsAsFactors = FALSE )
  
  if (!is.null(density.cutoff)){
    from <- as.numeric(unlist(strsplit(density.cutoff,":")))[1]
    to <- as.numeric(unlist(strsplit(density.cutoff,":")))[2]
    mypats <- f1[which(f1$density >= from & f1$density <= to),]
  } else {
    cat("Applying quantile cutoff...")
    quant.cutoff <- as.numeric(quantile(f1$Freq, probs = c(0.96), na.rm=TRUE))
    mypats <- f1[which(f1$Freq >= quant.cutoff),]
  }
  mydf <- f2[f2$Pattern %in% mypats$Pattern,]
  fwrite(x=mydf, file=paste0(out.dir, "/vals_filtered.txt"), 
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}
