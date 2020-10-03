runMethimputeRegions <- function(samplefiles, Regionfiles, run.chr.separate=TRUE, context, out.dir) {
  
  df.obs <- list()
  df.sim <- list()
  merge.list <- vector(mode="list")
  # Read the sample file with filenames and file paths
  filelist <- fread(samplefiles, header=TRUE)
  
  if (run.chr.separate==TRUE){
    
    for (i in 1:length(filelist$file)){
      for (j in 1:length(context)){
        methfn <- gsub(".*methylome_|\\_All.txt$", "", filelist$file[i])
        Regfiles <- list.files(Regionfiles, pattern=paste0("_", context[j], ".Rdata"), full.names = TRUE)
        for (k1 in 1:length(Regfiles)){
          tmp <- gsub(".*_regions_|\\.Rdata$", "", Regfiles[k1])
          chr <- gsub(paste0("_", context[j]), "", tmp)
          cat(paste0("Running file: ", methfn, " for context: ", context[j], " and chromosome: ", chr,"\n"), sep = "")
          out <- makeMethimpute(
            df=filelist$file[i],
            context=context[j],
            refRegion=Regfiles[k1],
            fit.plot=FALSE,
            merge.chr=TRUE,
            include.intermediate=FALSE, 
            probability="constrained",
            out.dir=out.dir,
            fit.name=paste0(methfn, "_", context[j], "_", chr),
            name=methfn)
          merge.list[[k1]] <- data.frame(out)
        }
        full.list <- rbindlist(merge.list)
        cat(paste0("Writing output to file", "\n"), sep="")
        modifiedExportMethylome(model=full.list, out.dir=out.dir, context=context[j], name=methfn)
      } 
    }
  } else {
    for (i in 1:length(filelist$file)){
      for (j in 1:length(context)){
        methfn <- gsub(".*methylome_|\\_All.txt$", "", filelist$file[i])
        Regfiles <- list.files(Regionfiles, pattern=paste0("_", context[j], ".Rdata"), full.names = TRUE)
        cat(paste0("Merging individual chr data for file: ", methfn, " context ", context[j], " .........\n"), sep="")
        for (k2 in 1:length(Regfiles)){
          f.file <- dget(Regfiles[k2])
          if (NROW(f.file$reg.obs)==0) {
            cat(paste0("Empty file ", basename(Regfiles[k2]), " .........\n"), sep="")
          } else {
            f.file$reg.obs$chr <- unlist(f.file$reg.obs$chr)
            f.file$reg.sim$chr <- unlist(f.file$reg.sim$chr)
            df.obs[[k2]] <- as.data.frame(f.file$reg.obs)
            df.sim[[k2]] <- as.data.frame(f.file$reg.sim)
          }
        }
        outlist <- list(reg.obs=do.call(rbind,df.obs),
                        reg.sim=do.call(rbind,df.sim),
                        context=context[j])
        
        regMerged <- paste(out.dir, "/Regions_merged_", context[j], ".Rdata", sep="")
        dput(outlist, regMerged)
        
        makeMethimpute(
          df=filelist$file[i],
          context=context[j],
          refRegion=regMerged,
          fit.plot=FALSE,
          merge.chr=FALSE,
          include.intermediate=FALSE, 
          probability="constrained",
          out.dir=out.dir,
          fit.name=paste0(methfn, "_", context[j]),
          name=methfn)
      }
    }
  }
}

#----------------------------------------------------------------
# mergeRegions <- function(samplefiles, context, input.dir, out.dir) {
#   merge.list <- list()
#   # Read the sample file with filenames
#   samplelist <- fread(samplefiles, header=T)
#   for (i in 1:length(samplelist)){
#     for (j in  1:length(context)){
#       mynames <- gsub(".*methylome_|\\_All.txt$", "", samplelist$file[i])
#       cat(paste0("Extracting files for ", mynames, " ", context[j], "\n"), sep = "")
#       mylist <- list.files(input.dir, pattern=mynames, full.names = TRUE)
#       if (length(mylist)!=0) {
#         selectfiles <- grep(context[j], mylist, value=TRUE)
#         print (selectfiles)
#         for (k in seq_along(selectfiles)){
#           f <- fread(selectfiles[k])
#           merge.list[[k]] <- f
#         }
#         full.list <- rbindlist(merge.list)
#         fwrite(x=full.list, file=paste0(out.dir, "/", mynames, "_", context[j],"_full.txt"), 
#                quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
#         rm(full.list)
#       } else {
#       cat(paste0("empty list ", "\n"), sep = "")
#       }
#     }
#   }
# }
