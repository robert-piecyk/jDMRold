#'
#' @param status.collect
#' @param rc.methlevel.collect
#' @param replicate.consensus
#' @importFrom data.table fread
#' @importFrom data.table rbindlist
#' @importFrom dplyr semi_join
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'


filterReplicateConsensus <- function(status.collect, rc.methlevel.collect, replicate.consensus, diff.ct=0.5){

  if (!is.null(status.collect$epiMAF)){
    status.collect <- status.collect[,-c("epiMAF")]
  }
  # deducing replicates info
  mycol <- names(status.collect)[4:NCOL(status.collect)]
  sampleinfo <- data.frame(do.call(rbind, strsplit(as.character(mycol),"_")))

  if (length(sampleinfo)==2){
    colnames(sampleinfo) <- c("sample", "replicate")
    dt <- data.frame()
    pb1 <- txtProgressBar(min = 1, max = NROW(status.collect), char = "=", style = 3, file = "")

    q <- lapply(1:NROW(status.collect), function(x){
      mypattern <- unlist(status.collect[x, 4:NCOL(status.collect)])
      df.bind <- cbind(sampleinfo, mypattern)
      for (m in unique(df.bind$sample)){
        total.reps <- length(df.bind$mypattern[df.bind$sample==m])
        rval <- round(replicate.consensus * total.reps)
        pattern.vals <- df.bind$mypattern[df.bind$sample==m]
        df.bind$diff.count[df.bind$sample==m] <- length(which(pattern.vals==1))/total.reps
        tt <- table(pattern.vals)
        if (max(tt) >= rval){
          df.bind$count[df.bind$sample==m] <- 0
        } else {
          df.bind$count[df.bind$sample==m] <- 1
        }
      }
      Sys.sleep(1/NROW(status.collect))
      setTxtProgressBar(pb1, x)
      close(pb1)

      df.gp <- group_by(df.bind, sample) %>% summarize(n = mean(diff.count))
      cb <- combn(df.gp$n,2)
      my.diff <- unlist(lapply(cb, function(x) abs(cb[1]-cb[2])))

      # allowing 50% difference between control and treatment groups
      if ((min(my.diff) >= diff.ct) && (sum(df.bind$count)==0)) {
        dt <- rbind(dt, status.collect[x,])
      }
      #print(df.bind)
      #if (sum(df.bind$count)==0) {
      #  dt <- rbind(dt, status.collect[x,])
      #}
    })
    out <- q[!sapply(q,is.null)]
    #status.collect <- q[-which(sapply(q, is.null))]
    df.status.collect <- rbindlist(out)
    if (NROW(df.status.collect) !=0){
      df.rc.methlevel.collect <- rc.methlevel.collect %>% dplyr::semi_join(df.status.collect,
                                                                           by=c("seqnames","start","end"))
      return(list(df.status.collect, df.rc.methlevel.collect))
    } else {
      message("\nEmpty dataframe. Nothing to write!")
      return(NULL)
    }
  } else {
    stop("Column for replicates is missing!!!")
  }
}

#------------------------------------------------------------------------------------------------

#' @param mat1
#' @param mat2
#' @param epiMAF
#' @importFrom dplyr semi_join
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'

filterEpiMAF <- function(mat1, mat2, epiMAF){

  pb2 <- txtProgressBar(min = 1, max = NROW(mat1), char = "=", style = 3, file = "")
  mat1$epiMAF <- 0

  for (i1 in 1:NROW(mat1)){
    mypattern <- unlist(mat1[i1, 4:(NCOL(mat1)-1)])
    mycount <- table(mypattern)
    epiMAF.out <- min(mycount)/length(mypattern)
    mat1$epiMAF[i1] <- floorDec(as.numeric(as.character(epiMAF.out)),5)

    Sys.sleep(1/NROW(mat1))
    setTxtProgressBar(pb2, i1)
  }
  close(pb2)

  df.status.collect <- mat1[which(mat1$epiMAF < epiMAF),]
  if (NROW(df.status.collect) !=0){
    df.rc.methlevel.collect <- mat2 %>% dplyr::semi_join(df.status.collect, by=c("seqnames","start","end"))
    return(list(df.status.collect, df.rc.methlevel.collect))
  } else {
    message("\nEmpty dataframe. Nothing to write!")
    return(NULL)
  }
}

#------------------------------------------------------------------------------------------------
#'
#' @param rcmethlvl
#' @param statecalls
#' @param gap
#' @import GenomicRanges
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
merge.bins <- function(rcmethlvl, statecalls, rc.methlvl.out){
  gap=1

  mylist <- list()
  result <- list()

  matrix1 <- as.data.frame(statecalls)
  matrix2 <- as.data.frame(rcmethlvl)

  # extract unique pattern
  extract.pattern <- unique(matrix1[,(4:NCOL(matrix1))])
  extract.pattern$pattern <- apply(extract.pattern[,c(1:NCOL(extract.pattern))], 1, paste, collapse="")

  # the state call matrix
  gr1 <- GRanges(seqnames=matrix1$seqnames, ranges=IRanges(start=matrix1$start, end=matrix1$end))
  values(gr1) <- cbind(values(gr1), pattern=apply(matrix1[,c(4:NCOL(matrix1))], 1, paste, collapse=""))

  message("\nNow, Merging overlapping and consecutive bins...\n")

  # this is for the state-calls: collapse overlapping bins if pattern is same
  grl_reduce <- unlist(GenomicRanges::reduce(split(gr1, gr1$pattern)))
  result <- sort(grl_reduce)
  result$pattern <- names(result)
  result <- data.frame(result)
  final.status.collect <- result %>% dplyr::left_join(extract.pattern, by=c("pattern"))

  if (rc.methlvl.out==TRUE){
    # the rcmethlvl matrix, also add the state-call pattern
    gr2 <- GRanges(seqnames=matrix2$seqnames, ranges=IRanges(start=matrix2$start, end=matrix2$end))
    values(gr2) <- cbind(values(gr2), values(gr1), DataFrame(matrix2[,c(4:NCOL(matrix2))]))

    # this is for the rcmethlvl: collapse bins and take average of the bins
    mycols <- colnames(matrix1)[4:NCOL(matrix1)]
    fn = function(u){
      out = GenomicRanges::reduce(u)
      for (x in 1:length(mycols)){
        eval(parse(text=paste0("out$", mycols[x], " = mean(u$", mycols[x], ")")))
      }
      return(out)
    }
    # split rcmthlvl matrix based on pattterns
    grl <- split(gr2, gr2$pattern)

    pb3 <- txtProgressBar(min = 1, max = length(grl), char = "=", style = 3, file = "")

    for (x in 1:length(grl)){
      a <- data.frame(grl[[x]])
      a1 <- a %>% arrange(pattern, start) %>% group_by(pattern) %>% mutate(indx = cumsum(start > lag(end, default = start[1]) + gap))
      a1.gr <- makeGRangesFromDataFrame(a1, keep.extra.columns=TRUE)
      df <- lapply(split(a1.gr, a1.gr$indx), fn)
      mylist[[x]] <- df

      Sys.sleep(0.05)
      setTxtProgressBar(pb3, x)
    }
    close(pb3)

    f.df <- unlist(mylist)
    f.df <- do.call(rbind, lapply(f.df, data.frame))
    f.df <- f.df[order(f.df[,1], f.df[,2]),]
    f.df[,(6:NCOL(f.df))] <- lapply(f.df[,(6:NCOL(f.df))], function(xy){ floorDec(xy,5) })
    final.rcmethlvl.collect <- subset(f.df, select = -c(strand))
    final.status.collect <- subset(final.status.collect, select = -c(strand, pattern))
    return(list(final.status.collect, final.rcmethlvl.collect))
  } else {
    final.status.collect <- subset(final.status.collect, select = -c(strand, pattern))
    final.rcmethlvl.collect <- NULL
    return(list(final.status.collect, final.rcmethlvl.collect))
  }
}
#------------------------------------------------------------------------------------------------

export.out <- function(out.rcmethlvl, out.statecalls, context, out.name1, out.name2, data.out){
  fwrite(x=out.statecalls,
         file=paste0(data.out, "/", context, "_", out.name1, ".txt"),
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  if (!is.null(out.rcmethlvl)) {
    fwrite(x=out.rcmethlvl,
           file=paste0(data.out, "/", context, "_", out.name2, ".txt"),
           quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  } else {
    message("Generate filtered matrix of recalibrated methylation levels set to FALSE!")
    message("------------------------------------------------------", "\n")
  }
}

#------------------------------------------------------------------------------------------------
#' filter DMR matrix
#'
#' Filters non-polymorphic patterns by default.
#' @param epiMAF.cutoff Filter for Minor Epi-Allele frequency (for e.g 0.33). Applicable for population level data. By default this option is set to NULL.
#' @param replicate.consensus Percentage of concordance in methylation states among samples with multiple replicates. Applicable for Control/treatment data. By default this option is set to NULL.
#' @param gridDMR Logical specifying if grid DMR approach was used for calling DMRs. By default this option is set to TRUE. If Cluster approach was used, set it to FALSE.
#' @param data.dir Directory containing DMR matrix files. Looks for files with suffix("_StateCalls.txt" and "_rcMethlvl.txt")
#' @param rc.methlvl.out Logical whether to output filtered matrix with methylation levels. By default this option is set to FALSE.
#' @param context.specific.DMRs Logical whether to output list of CG-only, CHG-only, CHH-only, nonCG and multi-context DMRs. By default this option is set to TRUE.
#' @importFrom data.table fread
#' @export
#'
filterDMRmatrix <- function(epiMAF.cutoff=NULL,
                            replicate.consensus=NULL,
                            gridDMR=TRUE,
                            data.dir,
                            rc.methlvl.out=FALSE,
                            samplefiles,
                            contexts=c("CG","CHG","CHH")) {

  ft <- fread(samplefiles)
  if (!is.null(ft$group)){
    ft$name <- paste0(ft$sample,"_", ft$replicate)
    gps <- ft$group[!ft$group %in% c('control')]
    gps <- unique(gps)
    list.collect1 <- list()
    list.collect2 <- list()
    for (m in seq_along(gps)){
      myvec <- c("control", gps[m])
      gp1 <- ft$name[which(ft$group==myvec[1])]
      gp2 <- ft$name[which(ft$group==myvec[2])]
      gp1.sample <- unique(ft$sample[which(ft$name==gp1)])
      gp2.sample <- unique(ft$sample[which(ft$name==gp2)])
      out.name <- paste0(gp1.sample, "_", gp2.sample)
      for (cn in seq_along(contexts)){
        fn1 <- paste0(data.dir, contexts[cn],"_", out.name ,"_StateCalls.txt")
        list.collect1[[cn]] <- fn1
      }
      list.collect2[[m]] <- list.collect1
    }
    list.status <- unlist(list.collect2)
  } else {
    list.status <- list.files(data.dir, pattern="_StateCalls.txt", full.names=TRUE)
  }
  if (length(list.status) != 0){
    for (i in seq_along(list.status)){
      context <- gsub("_StateCalls.txt", "", basename(list.status[i]))
      message("\nFiltering DMR matrix for ", context)

      #----------------------------------------------
      # Removing non-polymorphic/unchanged patterns
      #----------------------------------------------
      if (file.exists(list.status[i])){
        status.collect  <- fread(list.status[i], header=T)
      } else {
        stop("Files donot exist or is non-readable!")
      }

      message("Removing non-polymorphic patterns...")

      index <- which(rowSums(status.collect[,4:NCOL(status.collect)]) != 0 &
                       rowSums(status.collect[,4:NCOL(status.collect)]) != NCOL(status.collect)-3)

      status.collect <- status.collect[index,]

      rc.methlvl.name <- paste0(data.dir, "/", context, "_rcMethlvl.txt")
      rc.methlevel.collect <- fread(rc.methlvl.name, header=T)
      rc.methlevel.collect <- rc.methlevel.collect[index,]

      if (gridDMR==TRUE) {

        message("grid DMR set to TRUE")
        if (is.null(epiMAF.cutoff) && is.null(replicate.consensus)) {
          message("Both, epiMAF and replicate consensus set to NULL")
          out1=status.collect
          out2=rc.methlevel.collect
          out <- merge.bins(statecalls=out1, rcmethlvl=out2, rc.methlvl.out)
          export.out(out.statecalls=out[[1]],
                     out.rcmethlvl=out[[2]],
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        }
        #----------------------------------------------
        # Optional. Filtering out regions with epiMAF < Minor Allele Frequency
        #----------------------------------------------
        if (!is.null(epiMAF.cutoff)) {
          message("Filtering for epiMAF: ", epiMAF.cutoff, "\n")
          mydf <- filterEpiMAF(mat1=status.collect, mat2=rc.methlevel.collect, epiMAF=epiMAF.cutoff)
          if (!is.null(mydf)){
            # For Population data remove the epiMAF column
            out1=mydf[[1]][,-c("epiMAF")]
            out2=mydf[[2]]
            out <- merge.bins(statecalls=out1, rcmethlvl=out2, rc.methlvl.out)
            export.out(out.statecalls=out[[1]],
                       out.rcmethlvl=out[[2]],
                       context=context,
                       out.name1="StateCalls-filtered",
                       out.name2="rcMethlvl-filtered",
                       data.out=data.dir)
          }
        }
        #----------------------------------------------
        # Optional. Retaining samples based on replicate.consensus
        #----------------------------------------------
        if (!is.null(replicate.consensus)) {
          message("Filtering for replicate consensus...\n")
          mydf <- filterReplicateConsensus(status.collect, rc.methlevel.collect, replicate.consensus)
          if (!is.null(mydf)){
            out1=mydf[[1]]
            out2=mydf[[2]]
            out <- merge.bins(statecalls=out1, rcmethlvl=out2, rc.methlvl.out)
            export.out(out.statecalls=out[[1]],
                       out.rcmethlvl=out[[2]],
                       context=context,
                       out.name1="StateCalls-filtered",
                       out.name2="rcMethlvl-filtered",
                       data.out=data.dir)
          }
        }
      } else if (gridDMR==FALSE){

        message("\ngrid DMR set to FALSE. Running for Region DMRs")
        if (is.null(epiMAF.cutoff) && is.null(replicate.consensus)) {
          message("Both, epiMAF and replicate consensus set to NULL")
          out1=status.collect
          if (rc.methlvl.out==TRUE){
            out2=rc.methlevel.collect
          } else {
            out2 <- NULL
          }
          export.out(out.statecalls=out1,
                     out.rcmethlvl=out2,
                     context=context,
                     out.name1="StateCalls-filtered",
                     out.name2="rcMethlvl-filtered",
                     data.out=data.dir)
        }
        #----------------------------------------------
        # Optional. Filtering out regions with epiMAF < Minor Allele Frequency
        #----------------------------------------------
        if (!is.null(epiMAF.cutoff)) {
          message("Filtering for epiMAF: ", epiMAF.cutoff, "\n")
          mydf <- filterEpiMAF(mat1=status.collect, mat2=rc.methlevel.collect, epiMAF=epiMAF.cutoff)
          if (!is.null(mydf)){
            # For Population data remove the epiMAF column
            out1=mydf[[1]][,-c("epiMAF")]
            if (rc.methlvl.out==TRUE){
              out2=mydf[[2]]
            } else {
              out2 <- NULL
            }
            export.out(out.statecalls=out1,
                       out.rcmethlvl=out2,
                       context=context,
                       out.name1="StateCalls-filtered",
                       out.name2="rcMethlvl-filtered",
                       data.out=data.dir)
          }
        }
        #----------------------------------------------
        # Optional. Retaining samples based on replicate.consensus
        #----------------------------------------------
        if (!is.null(replicate.consensus)) {
          message("Filtering for replicate consensus...\n")
          mydf <- filterReplicateConsensus(status.collect, rc.methlevel.collect, replicate.consensus)
          if (!is.null(mydf)){
            out1=mydf[[1]]
            if (rc.methlvl.out==TRUE){
              out2=mydf[[2]]
            } else {
              out2 <- NULL
            }
            export.out(out.statecalls=out1,
                       out.rcmethlvl=out2,
                       context=context,
                       out.name1="StateCalls-filtered",
                       out.name2="rcMethlvl-filtered",
                       data.out=data.dir)
          }
        }
      }
    }
  } else {
    message("\nDMR matrix files do not exist!")
  }
  #context.specific.DMRs(data.dir)
}

#------------------------------------------------------------------------------------------------

DMR.list.out <- function(context.df, out.name, data.out){
  fwrite(x=context.df, file=paste0(data.out, "/", out.name,".txt"),
         quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}

#------------------------------------------------------------------------------------------------
#' @param data.dir
#' @importFrom IRanges subsetByOverlaps
#' @importFrom dplyr semi_join
#' @importFrom GenomicRanges findOverlaps
#' @importFrom data.table rbindlist
#' @importFrom data.table fread
#' @importFrom GenomicRanges intersect
#' @importFrom tidyr separate
#' @importFrom tidyr unnest
#' @importFrom dplyr mutate
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#'

extract.context.DMRs <- function(file1, file2, file3, tmp.name, data.dir){
  if (file.exists(file1) && file.exists(file2) && file.exists(file3)){
    CG.out <- data.table::fread(file1, header=TRUE)
    CHG.out <- data.table::fread(file2, header=TRUE)
    CHH.out <- data.table::fread(file3, header=TRUE)

    CG.gr <- GRanges(seqnames=CG.out$seqnames, ranges=IRanges(CG.out$start, end=CG.out$end))
    CHG.gr <- GRanges(seqnames=CHG.out$seqnames, ranges=IRanges(CHG.out$start, end=CHG.out$end))
    CHH.gr <- GRanges(seqnames=CHH.out$seqnames, ranges=IRanges(CHH.out$start, end=CHH.out$end))

    #CG.only
    message("Generating CG-only DMRs")
    out1 <- IRanges::subsetByOverlaps(CG.gr, CHG.gr, invert = TRUE)
    CG.1 <- IRanges::subsetByOverlaps(out1, CHH.gr, invert = TRUE)
    CG.1 <- as.data.frame(CG.1)
    CG.1$seqnames <- as.integer(as.character(CG.1$seqnames))
    CG.only <- CG.out %>% dplyr::semi_join(CG.1, by = c("seqnames","start","end"))
    DMR.list.out(context.df=CG.only,
                 out.name=paste0(tmp.name,"CG-only-DMRs"),
                 data.out=data.dir)
    message("Done!")

    #CHG.only
    message("Generating CHG-only DMRs")
    out2 <- IRanges::subsetByOverlaps(CHG.gr, CG.gr, invert = TRUE)
    CHG.1 <- IRanges::subsetByOverlaps(out2, CHH.gr, invert = TRUE)
    CHG.1 <- as.data.frame(CHG.1)
    CHG.1$seqnames <-as.integer(as.character(CHG.1$seqnames))
    CHG.only <- CHG.out %>% dplyr::semi_join(CHG.1, by = c("seqnames","start","end"))
    DMR.list.out(context.df=CHG.only,
                 out.name=paste0(tmp.name,"CHG-only-DMRs"),
                 data.out=data.dir)
    message("Done!")

    #CHH.only
    message("Generating CHH-only DMRs")
    out3 <- IRanges::subsetByOverlaps(CHH.gr, CG.gr, invert = TRUE)
    CHH.1 <- IRanges::subsetByOverlaps(out3, CHG.gr, invert = TRUE)
    CHH.1 <- as.data.frame(CHH.1)
    CHH.1$seqnames <-as.integer(as.character(CHH.1$seqnames))
    CHH.only <- CHH.out %>% dplyr::semi_join(CHH.1, by = c("seqnames","start","end"))
    DMR.list.out(context.df=CHH.only,
                 out.name=paste0(tmp.name, "CHH-only-DMRs"),
                 data.out=data.dir)
    message("Done!")

    #non-CG
    message("Generating non-CG DMRs")
    nonCG.collect <- list()
    overlaps.nonCG <- GenomicRanges::findOverlaps(CHG.gr, CHH.gr, ignore.strand=TRUE)
    overlaps.hits.nonCG <- IRanges::subsetByOverlaps(CHG.gr, CHH.gr)
    mcols(overlaps.hits.nonCG)$DMRs.CHH.coord<- CharacterList(split(CHH.gr[subjectHits(overlaps.nonCG)], queryHits(overlaps.nonCG)))
    out.nonCG <- IRanges::subsetByOverlaps(overlaps.hits.nonCG, CG.gr, invert = TRUE)
    nonCG <- data.frame(out.nonCG) %>% dplyr::mutate(DMRs.CHH.coord = strsplit(as.character(DMRs.CHH.coord), ",")) %>% tidyr::unnest(c(DMRs.CHH.coord))
    nonCG.clean <- data.frame(lapply(nonCG, function(k) gsub ("[\\c]|[()]|\"|^ .", "", k)))
    nonCG.clean <- nonCG.clean %>% tidyr::separate(DMRs.CHH.coord, c("CHH.seqnames","CHH.start","CHH.stop"), sep = '([-:])')
    nonCG.clean <- nonCG.clean[-c(4,5)]
    colnames(nonCG.clean)[1] <- "CHG.seqnames"
    colnames(nonCG.clean)[2] <- "CHG.start"
    colnames(nonCG.clean)[3] <- "CHG.stop"

    pb4 <- txtProgressBar(min = 1, max = NROW(nonCG.clean), char = "=", style = 3, file = "")

    for (i1 in 1:NROW(nonCG.clean)){
      myrow <- nonCG.clean[i1,]
      a=makeGRangesFromDataFrame(myrow[,c("CHG.seqnames","CHG.start","CHG.stop")])
      b=makeGRangesFromDataFrame(myrow[,c("CHH.seqnames","CHH.start","CHH.stop")])
      out <- data.frame(GenomicRanges::intersect(a,b))
      myrow$merged.seqnames <- out$seqnames
      myrow$merged.start <- out$start
      myrow$merged.stop <- out$end
      nonCG.collect[[i1]] <- data.frame(myrow)
      Sys.sleep(1/NROW(nonCG.clean))
      setTxtProgressBar(pb4, i1)
    }
    close(pb4)
    DMR.list.out(context.df=data.table::rbindlist(nonCG.collect),
                 out.name=paste0(tmp.name,"nonCG-DMRs"),
                 data.out=data.dir)
    message("Done!")

    #multi-context
    message("Generating multi-context DMRs")
    multi.context.collect <- list()
    overlaps.1 <- GenomicRanges::findOverlaps(CG.gr, CHG.gr, ignore.strand=TRUE)
    overlaps.hits.1 <- IRanges::subsetByOverlaps(CG.gr, CHG.gr)
    mcols(overlaps.hits.1)$DMRs.CHG.coord <- CharacterList(split(CHG.gr[subjectHits(overlaps.1)], queryHits(overlaps.1)))
    overlaps.2 <- GenomicRanges::findOverlaps(overlaps.hits.1, CHH.gr, ignore.strand=TRUE)
    overlaps.hits.2 <- IRanges::subsetByOverlaps(overlaps.hits.1, CHH.gr)
    if (NROW(overlaps.hits.2)!=0){
      mcols(overlaps.hits.2)$DMRs.CHH.coord <- CharacterList(split(CHH.gr[subjectHits(overlaps.2)], queryHits(overlaps.2)))
      multi.context.1 <- data.frame(overlaps.hits.2) %>% dplyr::mutate(DMRs.CHG.coord = strsplit(as.character(DMRs.CHG.coord), ",")) %>% tidyr::unnest(c(DMRs.CHG.coord))
      multi.context.1.clean <- data.frame(lapply(multi.context.1, function(k) gsub ("[\\c]|[()]|\"|^ .", "", k)))
      multi.context.1.clean <- multi.context.1.clean %>% tidyr::separate(DMRs.CHG.coord, c("CHG.seqnames","CHG.start","CHG.stop"), sep = '([-:])')
      multi.context.2.clean <- data.frame(multi.context.1.clean) %>% dplyr::mutate(DMRs.CHH.coord = strsplit(as.character(DMRs.CHH.coord), ",")) %>% tidyr::unnest(c(DMRs.CHH.coord))
      multi.context.2.clean <- multi.context.2.clean %>% tidyr::separate(DMRs.CHH.coord, c("CHH.seqnames","CHH.start","CHH.stop"), sep = '([-:])')
      multi.context.2.clean <- multi.context.2.clean[-c(4,5)]
      colnames(multi.context.2.clean)[1] <- "CG.seqnames"
      colnames(multi.context.2.clean)[2] <- "CG.start"
      colnames(multi.context.2.clean)[3] <- "CG.stop"

      pb5 <- txtProgressBar(min = 1, max = NROW(multi.context.2.clean), char = "=", style = 3, file = "")

      for (i2 in 1:NROW(multi.context.2.clean)){
        myrow.x <- multi.context.2.clean[i2,]
        a=makeGRangesFromDataFrame(myrow.x[,c("CG.seqnames","CG.start","CG.stop")])
        b=makeGRangesFromDataFrame(myrow.x[,c("CHG.seqnames","CHG.start","CHG.stop")])
        c=makeGRangesFromDataFrame(myrow.x[,c("CHH.seqnames","CHH.start","CHH.stop")])
        out1 <- GenomicRanges::intersect(a,b)
        out2 <- data.frame(GenomicRanges::intersect(out1,c))
        if(NROW(out2)!=0){
          myrow.x$merged.seqnames <- out2$seqnames
          myrow.x$merged.start <- out2$start
          myrow.x$merged.stop <- out2$end
          multi.context.collect[[i2]] <- data.frame(myrow.x)
        }
        Sys.sleep(1/NROW(multi.context.2.clean))
        setTxtProgressBar(pb5, i2)
      }
      close(pb5)
      f <- data.table::rbindlist(multi.context.collect)

      DMR.list.out(context.df=f,
                   out.name=paste0(tmp.name, "multi-context-DMRs"),
                   data.out=data.dir)
    } else {
      message("No multi-context DMRs found!")
    }
    message("Done!")
  } else {
    stop("Filtered DMR matrix files for all contexts donot exist!")
  }
}

context.specific.DMRs <- function(samplefiles, data.dir){
  ft <- fread(samplefiles)
  if (!is.null(ft$group)){
    ft$name <- paste0(ft$sample,"_", ft$replicate)
    gps <- ft$group[!ft$group %in% c('control')]
    gps <- unique(gps)
    for (m in seq_along(gps)){
      myvec <- c("control", gps[m])
      gp1 <- ft$name[which(ft$group==myvec[1])]
      gp2 <- ft$name[which(ft$group==myvec[2])]

      gp1.sample <- unique(ft$sample[which(ft$name==gp1)])
      gp2.sample <- unique(ft$sample[which(ft$name==gp2)])
      message("Generating context specific DMRs for ", gp1.sample, "-", gp2.sample,"\n")
      CG.f <- paste0(data.dir, "CG_", gp1.sample, "_", gp2.sample, "_StateCalls-filtered.txt")
      CHG.f <- paste0(data.dir, "CHG_", gp1.sample, "_", gp2.sample, "_StateCalls-filtered.txt")
      CHH.f <- paste0(data.dir, "CHH_", gp1.sample, "_", gp2.sample, "_StateCalls-filtered.txt")
      extract.context.DMRs(file1=CG.f,
                           file2=CHG.f,
                           file3=CHH.f,
                           tmp.name=paste0(gp1.sample, "_", gp2.sample, "_"),
                           data.dir=data.dir)
    }
  } else {
    message("Generating context specific DMRs. No groups found!\n")
    output <- extract.context.DMRs(file1=paste0(data.dir,"CG_StateCalls-filtered.txt"),
                                   file2=paste0(data.dir,"CHG_StateCalls-filtered.txt"),
                                   file3=paste0(data.dir,"CHH_StateCalls-filtered.txt"),
                                   tmp.name="",
                                   data.dir=data.dir)
  }
}



