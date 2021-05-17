# This function will merge (column 6) state calls and (column 7) rc.meth.lvl from all samples into one dataframe
# makes list of 2 dataframes
merge.cols <- function(filepath, colm, include.intermediate) {

  mylist <- list()
  for (l in 1:length(colm)){
    extract <- lapply(filepath, function(k){
      f <- fread(k, header=FALSE, skip=1, select=c(1, 2, 3, colm[l]))
      if (colm[l]==6) {
        if (include.intermediate==TRUE) {
          f[,4] <- ifelse(f[,4] == "U", yes = 0, (ifelse(f[,4] == "I", yes = 0.5, no = 1)))
        } else { 
          f[,4] <- ifelse(f[,4] == "U", yes = 0, no = 1)
        }
      }
      colnames(f)[4] <- basename(k)
      return(f)
      
    })
    df <- Reduce(function(x, y) {
      dplyr::inner_join(x, y, by=c("V1","V2","V3"))
    }, extract)
    
    mylist[[l]] <- df
  }
  return(mylist)
}

makeDMRmatrix <- function(context, samplefiles, input.dir, out.dir, include.intermediate) {
  # Read the sample file with filenames
  samplelist <- fread(samplefiles, header=T)
  for (j in  1:length(context)){
    
    # list all files in the input directory
    extractflist <- list.files(input.dir, pattern=paste0(context[j],".txt"), full.names=TRUE)
    
    #extract file subsets for construction of DMRmatrix
    if (length(extractflist) != 0){
      mynames <- gsub(paste0("_", context[j], ".txt$"), "", basename(extractflist))
      selectlist <- list()
      message("\nExtracting filenames and matching them....\n")
      for (a1 in seq_along(mynames)){
        as <- samplelist[grepl(paste0("_",mynames[a1],"_"), samplelist$file),]
        if (NROW(as)==1){
          as$full.path.MethReg <- grep(paste0("/", mynames[a1], "_", context[j], ".txt", sep=""), extractflist, value=TRUE)
          message("\nRegion file ", basename(as$full.path.MethReg)," found !")
          selectlist[[a1]] <- as
        } else { 
          message("\nMultiple files with string match ", mynames[a1]," found !")
        }
      }
      flist <- rbindlist(selectlist)
      #print(flist)
      
      # Assign unique names for samples with or without replicate data
      if (!is.null(flist$replicate)) {
        message(paste0("\nRunning context ", context[j], ". Input data with replicates, creating unique sample names...\n"), sep = "")
        flist$name <- paste0(flist$sample,"_", flist$replicate)
      } else {
        flist$name <- flist$sample 
      }
      
      message(paste0("Now, constructing DMR matrix for ", context[j]), sep = "")
      
      # merge samples by Chr coordinates
      #(column 6) state-calls and (column 7) rc.meth.lvl
      mydf <- merge.cols(filepath=flist$full.path.MethReg, include.intermediate=include.intermediate, colm=c(5, 6, 7))
      
     # list containing state calls
      status.collect <- mydf[[2]]
      # renaming file names with sample names
      for (a in 4:length(colnames(status.collect))) {
        for (n in 1:length(flist$name)) {
          if (colnames(status.collect)[a] == basename(flist$full.path.MethReg)[n]) {
            colnames(status.collect)[a] = flist$name[n]
          }
        }
      }
      # list containing rcmethlvls
      rc.methlevel.collect <- mydf[[3]]
      # renaming file names with sample names
      for (a in 4:length(colnames(rc.methlevel.collect))) {
        for (n in 1:length(flist$name)) {
          if (colnames(rc.methlevel.collect)[a] == basename(flist$full.path.MethReg)[n]) {
            colnames(rc.methlevel.collect)[a] = flist$name[n]
          }
        }
      }
      
       # list containing postmax
      postMax.collect <- mydf[[1]]
      # renaming file names with sample names
      for (a in 4:length(colnames(postMax.collect))) {
        for (n in 1:length(flist$name)) {
          if (colnames(postMax.collect)[a] == basename(flist$full.path.MethReg)[n]) {
            colnames(postMax.collect)[a] = flist$name[n]
          }
        }
      }

      names(status.collect)[1] <- "seqnames"
      names(status.collect)[2] <- "start"
      names(status.collect)[3] <- "end"
      
      names(rc.methlevel.collect)[1] <- "seqnames"
      names(rc.methlevel.collect)[2] <- "start"
      names(rc.methlevel.collect)[3] <- "end"

      names(postMax.collect)[1] <- "seqnames"
      names(postMax.collect)[2] <- "start"
      names(postMax.collect)[3] <- "end"
      
      fwrite(x=status.collect, file=paste0(out.dir,"/", context[j],"_StateCalls.txt"),
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=rc.methlevel.collect, file=paste0(out.dir,"/", context[j],"_rcMethlvl.txt"), 
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
      fwrite(x=postMax.collect, file=paste0(out.dir,"/", context[j],"_postMax.txt"), 
             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
      message(paste0("\nDone! \n"), sep = "")
    } else{
      message(paste0("Files for context ",context[j]," do not exist\n"), sep="")
    }
  }
}