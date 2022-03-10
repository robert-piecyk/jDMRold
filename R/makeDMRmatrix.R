
#' @param filepath
#' @param colm
#' @param include.intermediate
#' @param mincov
#' @param nCytosines
#' @importFrom dplyr inner_join
#' @export
#'

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

write.out <- function(out.df, data.dir, out.name, contexts){
  fwrite(x=out.df, file=paste0(data.dir,"/", contexts,"_",out.name,".txt"),
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

#' Builds a DMR matrix for all samples
#'
#' This function generates a binary matrix, a matrix of recalibrated methylation levels and posterior probabilities for all samples.
#'
#' @param context cytosine context
#' @param samplefiles file containing full path of base level methylation calls, sample names and replicates(optional)
#' @param include.intermediate A logical specifying whether or not the intermediate component should be included in the HMM model.By default this option is set as FALSE.
#' @param input.dir input directory containing all region level methylome calls
#' @param out.dir output directory
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table rbindlist
#' @export
#'
makeDMRmatrix <- function(contexts=c("CG","CHG","CHH"), postMax.out=FALSE, samplefiles, input.dir, out.dir, include.intermediate=FALSE) {
  # Read the sample file with filenames
  samplelist <- fread(samplefiles, header=T)
  for (j in  1:length(contexts)){

    # list all files in the input directory
    extractflist <- list.files(input.dir, pattern=paste0(contexts[j],".txt"), full.names=TRUE)

    #extract file subsets for construction of DMRmatrix
    if (length(extractflist) != 0){
      mynames <- gsub(paste0("_", contexts[j], ".txt$"), "", basename(extractflist))
      selectlist <- list()
      message("\nExtracting filenames and matching them....")
      for (a1 in seq_along(mynames)){
        pat1 <- paste0("_",mynames[a1],"_","|","_",mynames[a1])
        #pat1 <- paste0("_",mynames[a1],"_")
        as <- samplelist[grepl(pat1, samplelist$file),]
        if (NROW(as)==1){
          as$full.path.MethReg <- grep(paste0("/", mynames[a1], "_", contexts[j], ".txt", sep=""), extractflist, value=TRUE)
          message("\n", basename(as$full.path.MethReg)," found !")
          selectlist[[a1]] <- as
        } else {
          message("\nMultiple files with string match ", mynames[a1]," found !")
        }
      }
      flist <- data.table::rbindlist(selectlist)
      #print(flist)

      # Assign unique names for samples with or without replicate data
      if (!is.null(flist$replicate)) {
        message(paste0("\nRunning context ", contexts[j], ". Input data with replicates, creating unique sample names...\n"), sep = "")
        flist$name <- paste0(flist$sample,"_", flist$replicate)
      } else {
        flist$name <- flist$sample
      }

      message(paste0("\nNow, constructing DMR matrix for ", contexts[j]), sep = "")

      # merge samples by Chr coordinates
      #(column 6) state-calls and (column 7) rc.meth.lvl
      mydf <- merge.cols(filepath=flist$full.path.MethReg,
                         include.intermediate=include.intermediate,
                         colm=c(5, 6, 7))

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
      if (postMax.out==TRUE){
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
        names(postMax.collect)[1] <- "seqnames"
        names(postMax.collect)[2] <- "start"
        names(postMax.collect)[3] <- "end"
        write.out(out.df=postMax.collect, data.dir=out.dir, out.name="postMax",contexts=contexts[j])
      }

      names(status.collect)[1] <- "seqnames"
      names(status.collect)[2] <- "start"
      names(status.collect)[3] <- "end"

      names(rc.methlevel.collect)[1] <- "seqnames"
      names(rc.methlevel.collect)[2] <- "start"
      names(rc.methlevel.collect)[3] <- "end"


      write.out(out.df=status.collect, data.dir=out.dir, out.name="StateCalls",contexts=contexts[j])
      write.out(out.df=rc.methlevel.collect, data.dir=out.dir, out.name="rcMethlvl",contexts=contexts[j])


      message(paste0("\nDone! \n"), sep = "")
    } else{
      message(paste0("Files for context ",contexts[j]," do not exist\n"), sep="")
    }
  }
}

split.groups <- function(samplefiles, postMax.out=FALSE, contexts=c("CG","CHG","CHH"), input.dir, out.dir){
  ft <- fread(samplefiles)
  ft$name <- paste0(ft$sample,"_", ft$replicate)
  gps <- ft$group[!ft$group %in% c('control')]
  gps <- unique(gps)
  for (m in seq_along(gps)){
    myvec <- c("control", gps[m])
    gp1 <- ft$name[which(ft$group==myvec[1])]
    gp2 <- ft$name[which(ft$group==myvec[2])]
    gp1.sample <- unique(ft$sample[which(ft$name==gp1)])
    gp2.sample <- unique(ft$sample[which(ft$name==gp2)])
    out.name <- paste0(gp1.sample, "_", gp2.sample)
    for (cn in seq_along(contexts)){
      fn1 <- paste0(input.dir, contexts[cn], "_StateCalls.txt")
      if (file.exists(fn1)) {
        StateCall <- fread(fn1)
        df1 <- dplyr::select(StateCall, seqnames, start, end, gp1, gp2)
        fwrite(x=df1,
               file=paste0(out.dir, "/", contexts[cn], "_", out.name, "_StateCalls.txt"),
               quote=FALSE,
               row.names=FALSE,
               col.names=TRUE,
               sep="\t")
      }

      fn2 <- paste0(input.dir, contexts[cn], "_rcMethlvl.txt")
      if (file.exists(fn2)) {
        rcMethlvl <- fread(fn2)
        df2 <- dplyr::select(rcMethlvl, seqnames, start, end, gp1, gp2)
        fwrite(x=df2,
               file=paste0(out.dir, "/", contexts[cn],"_", out.name, "_rcMethlvl.txt"),
               quote=FALSE,
               row.names=FALSE,
               col.names=TRUE,
               sep="\t")
      }
      if (postMax.out==TRUE){
        fn3 <- paste0(input.dir, contexts[cn], "_postMax.txt")
        if (file.exists(fn3)) {
          postMax <- fread(fn3)
          df3 <- dplyr::select(postMax, seqnames, start, end, gp1, gp2)
          fwrite(x=df3,
                 file=paste0(out.dir, "/", contexts[cn], "_", out.name, "_postMax.txt"),
                 quote=FALSE,
                 row.names=FALSE,
                 col.names=TRUE,
                 sep="\t")
        }
      }
    }
  }
}
