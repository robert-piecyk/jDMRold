
makeDMRmatrix <- function(context, chr, samplefile, out.dir) {
  
  mylist <- list()
  
  for (i in  1:length(context)){
    for (j in 1:length(chr)){
      cat(paste0("Running for ", chr[j], "_", context[i], "\n"), sep = "")
      
      # extract samples in listFiles.fn
      samplelist <- fread(samplefile, header=T)
      extractp <- paste0(chr[j],"_",context[i])
      samplelist <- samplelist[grep(extractp, samplelist$file),]
      print (samplelist$file)
      
      # (column 6) state calls and (column 7) rc.meth.lvl for all samples and make 2 matrices
      merge.cols <- function(colm) {
        for (l in 1:length(colm)){
          extract <- lapply(samplelist$file, function(k){
            f <- fread(k, header=FALSE, skip=1, 
                       select=c(1, 2, 3, colm[l]))
            if (colm[l]==6) {
              f[,4] <- ifelse(f[,4] == "U", yes = 0, no = 1)
            }
            colnames(f)[4] <- basename(k)
            return(f)
          })
          df <- Reduce(function(x, y) {
            dplyr::inner_join(x, y, by=c("V1","V2","V3"))
          }, extract)
          
          # renaming file name with a sample name
          for (a in 4:length(colnames(df))) {
            for (n in 1:length(samplelist$sample)) {
              if (colnames(df)[a] == basename(samplelist$file)[n]) {
                colnames(df)[a] = samplelist$sample[n]
              }
            }
          }
          mylist[[l]] <- df
        }
        return(mylist)
      }
      
      #(column 6) state-calls and (column 7) rc.meth.lvl
      mydf <- merge.cols(colm=c(6, 7))
      
      # list containing the state calls
      status.collect <- mydf[[1]]
      # list containing the rcmethlvls
      rc.methlevel.collect <- mydf[[2]]
      
      # removing unchanged patterns
      index <- which(rowSums(status.collect[,(4:NCOL(status.collect))]) != 0 & 
                       rowSums(status.collect[,(4:NCOL(status.collect))]) != (NCOL(status.collect)-3))
      
      status.collect <- status.collect[index,]
      rc.methlevel.collect <- rc.methlevel.collect[index,]
      
      names(status.collect)[1] <- "seqnames"
      names(status.collect)[2] <- "start"
      names(status.collect)[3] <- "end"
      
      names(rc.methlevel.collect)[1] <- "seqnames"
      names(rc.methlevel.collect)[2] <- "start"
      names(rc.methlevel.collect)[3] <- "end"
      
      fwrite(x=status.collect, file=paste0(out.dir, chr[j], "_", context[i],"_state-calls.txt"),
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=rc.methlevel.collect, file=paste0(out.dir, chr[j], "_", context[i],"_rcmethlvl.txt"), 
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    }
  }
}
