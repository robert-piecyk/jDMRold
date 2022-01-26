#' @param model
#' @param out.dir
#' @param name
#' @param distcor
#' @param skip
#' @param plot.parameters
#' @param df
#' @param refRegion
#' @param context Cytosine context
#' @param fit.plot
#' @param fit.name
#' @param refRegion
#' @param include.intermediate
#' @param probability
#' @param mincov
#' @param nCytosines
#' @import ggplot2
#' @importFrom minpack.lm nlsLM
#' @importFrom grDevices dev.off pdf
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom stringr str_replace_all
#' @import dplyr
#' @importFrom stats na.omit coefficients
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom methimpute callMethylation
#' @importFrom methimpute distanceCorrelation
#' @export

modified.estimateTransDist <- function(distcor, skip=2, plot.parameters=TRUE) {

  ## Context correlation fits and plots
  contexts <- dimnames(distcor$data)[[1]]
  cor.array <- distcor$data
  maxweights <- numeric()
  params.list <- list()
  miny <- min(cor.array, na.rm = TRUE)
  dfs <- list()
  for (c1 in 1:length(contexts)) {
    for (c2 in 1:length(contexts)) {
      context.transition <- paste0(contexts[c1], '-', contexts[c2])
      if (distcor$separate.contexts) {
        if (c1 != c2) {
          next
        }
      }
      if (c1 <= c2) {
        df <- data.frame(distance = as.numeric(dimnames(cor.array)[[3]]),
          correlation = cor.array[c1,c2,,'correlation'],
          weight = cor.array[c1,c2,,'weight'],
          from = contexts[c1], to = contexts[c2])

        ## Fit
        y <- df$correlation[(skip+1):nrow(df)]
        x <- df$distance[(skip+1):nrow(df)]
        weight <- df$weight[(skip+1):nrow(df)]
        startvalues <- list(a0 = stats::na.omit(y)[1], D = 50)
        p <- NULL

        if (is.null(p)) {
          startvalues <- list(a0 = stats::na.omit(y)[1])
          p <- tryCatch({
            fit <- minpack.lm::nlsLM(y ~ a0 * exp(-x/Inf), start=startvalues, weights=weight)
            s <- summary(fit)
            c <- stats::coefficients(s)
            params <- c[1:length(startvalues)]
            names(params) <- names(startvalues)
            params <- as.list(params)
            params$D <- Inf
            params
          }, error = function(e) {
            startvalues$D <- Inf
            startvalues
          })
        }

        ## Check if we have negative D
        if (p$D <= 0) {
          p$D <- Inf
        }
        params.list[[context.transition]] <- p

        ## Plot
        df$correlation.fit <- p$a0 * exp(-df$distance/p$D)
        df$logweight <- log(df$weight+1)
        dfs[[context.transition]] <- df
        maxweights[context.transition] <- max(df$logweight, na.rm = TRUE)
      }
    }
  }
  maxweight <- max(maxweights, na.rm = TRUE)

  ## Plot correlation
  df <- do.call(rbind, dfs)
  df$a0 <- round(sapply(params.list[paste0(df$from, '-', df$to)], '[[', 'a0'), 2)
  df$D <- round(sapply(params.list[paste0(df$from, '-', df$to)], '[[', 'D'), 0)
  df$params <- paste0("a0 = ", df$a0, ", D = ", df$D)
  ggplt <- ggplot(df) + theme_bw() + geom_line(aes_string(x='distance', y='correlation', alpha='logweight'))
  ggplt <- ggplt + geom_line(aes_string(x='distance', y='correlation.fit'), col='blue')
  if (plot.parameters) {
    ggplt <- ggplt + geom_text(aes_string(label='params'), x=max(df$distance, na.rm = TRUE), y=max(df$correlation, na.rm = TRUE), vjust=1, hjust=1)
  }
  ggplt <- ggplt + xlab('distance in [bp]')
  ggplt <- ggplt + facet_grid(from ~ to)
  if (miny < 0) {
    ggplt <- ggplt + geom_hline(aes_string('yintercept'=0), linetype=2, alpha=0.5)
  }

  transDist <- sapply(params.list, '[[', 'D')
  return(list(transDist=transDist, plot=ggplt))
}

#--------------------------------------------------------------------------
modifiedExportMethylome <- function(model, out.dir, context, name) {
    #data <- model$data
    data <- model
    final_dataset <- as(data, 'data.frame')
    final_dataset <- final_dataset[,c('seqnames','start','end','strand',
      'context','counts.methylated','counts.total',
      'posteriorMax','status','rc.meth.lvl')]

    # dropping columns
    drops <- c('width','strand','clusterlen','counts.methylated',
      'counts.total', 'distance', 'transitionContext', 'posteriorMeth','posteriorUnmeth')
    final_dataset <- final_dataset[ , !(names(final_dataset) %in% drops)]
    #------------------------------------------------------------------
    # convert full string into M/U/I
    final_dataset <- statusStringCheck(final_dataset)
    #------------------------------------------------------------------
    # take 4 digit of decimal value posteriorMax column
    final_dataset$posteriorMax <-floorDec(as.numeric(as.character(final_dataset$posteriorMax)),5)
    final_dataset$rc.meth.lvl <- floorDec(as.numeric(as.character(final_dataset$rc.meth.lvl)),5)
    final_dataset$seqnames <- as.character(final_dataset$seqnames)

    saveFile <- paste0(out.dir, "/", name, "_", context, ".txt")
    fwrite(final_dataset, file = saveFile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
    return (final_dataset)
}

#--------------------------------------------------------------------------

makeRegionsImpute <- function(df, context, refRegion, mincov, nCytosines) {

  #regions file
  #tmp_reg <- dget(refRegion)
  tmp_reg <- refRegion
  data <- as.data.frame(tmp_reg$reg.obs)
  data <- data %>% dplyr::filter(data$chr != "M" & data$chr != "C")
  #colnames(data)[which(names(data) == "cluster.size")] <- "cluster.length"

  #reference methimpute file
  ref_data <- data.table::fread(df, skip = 1, select = c("V1","V2","V3","V4","V5","V6"))

  #remove Mt and chloroplast coordinates. Following is for Arabidopsis only
  ref_data <- ref_data %>% dplyr::filter(ref_data$V1 != "M" & ref_data$V1 != "C")

  #filtering by coverage
  if (mincov>0){
    cat(paste0("Filtering for coverage: ", mincov,"\n"), sep = "")
    ref_data <- ref_data[which(ref_data$V6 >= mincov),]
  }

  ref_data <- ref_data[which(ref_data$V4==context),]

  data_gr <- GRanges(seqnames=data$chr,
                     ranges=IRanges(start=data$start, end=data$end),
                     clusterlen=data$cluster.length,
                     context=as.factor(context))

  ref_gr <- GRanges(seqnames=ref_data$V1,
                    ranges=IRanges(start=ref_data$V2, width=1),
                    context=as.factor(context),
                    methylated=ref_data$V5,
                    total=ref_data$V6)

  data_gr$cytosineCount <- GenomicRanges::countOverlaps(data_gr, ref_gr)

  #filtering by number of cytosines
  if (nCytosines>0){
    cat(paste0("Minimum number of cytosines per region: ", nCytosines,"\n"), sep = "")
    data_gr <- data_gr[which(data_gr$cytosineCount >= nCytosines),]
  }

  counts <- array(NA, dim=c(length(data_gr), 2), dimnames=list(NULL, c("methylated", "total")))

  overlaps <- IRanges::findOverlaps(ref_gr, data_gr)

  overlaps.hits <- ref_gr[S4Vectors::queryHits(overlaps)]
  if (NROW(overlaps.hits) != 0){

    methylated <- stats::aggregate(overlaps.hits$methylated, list(S4Vectors::subjectHits(overlaps)), FUN=sum)
    total <- stats::aggregate(overlaps.hits$total, list(S4Vectors::subjectHits(overlaps)), FUN=sum)

    if (NROW(methylated) != NROW(counts) ){
      missingr <- which(!rownames(data.frame(data_gr)) %in% methylated$Group.1)

      methylated <- dplyr::bind_rows(data.frame(Group.1=missingr, x=0), methylated)
      methylated <- methylated[order(methylated$Group.1),]

      total <- dplyr::bind_rows(data.frame(Group.1=missingr,x=0), total)
      total <- total[order(total$Group.1),]

      # for(item in seq_len(NROW((missingr)))){
      #   methylated <- rbind(c(missingr[item],0), methylated)
      #   methylated <- methylated[order(methylated$Group.1),]
      #   total <- rbind(c(missingr[item],0), total)
      #   total <- total[order(total$Group.1),]
      # }
    }
    counts[,"methylated"] <- methylated$x
    counts[,"total"] <- total$x
    data_gr$counts <- counts
  }
  rm(ref_data, ref_gr, overlaps, overlaps.hits)
  return(data_gr)
}

#--------------------------------------------------------------------------
makeMethimpute <- function(df, context, fit.plot, fit.name, refRegion,
                         include.intermediate, probability, out.dir, name, mincov, nCytosines){
  methylome.data <- makeRegionsImpute(df, context, refRegion, mincov, nCytosines)
  if (!is.null(methylome.data$counts)) {
    quant.cutoff <- as.numeric(quantile(methylome.data$counts[,"total"], probs = c(0.96), na.rm=TRUE))
    distcor <- distanceCorrelation(methylome.data, distances=0:100)
    fit <- modified.estimateTransDist(distcor)

    if (fit.plot==TRUE){
      print(paste0("Generating fit plot...", name))
      pdf(paste0(out.dir, "/", fit.name, "-fit.pdf", sep = ""))
      print(fit)
      dev.off()
    }

    model <- callMethylation(data = methylome.data,
                             transDist = fit$transDist,
                             count.cutoff = quant.cutoff,
                             max.time = Inf,
                             max.iter = Inf,
                             include.intermediate = include.intermediate,
                             update = probability)

    methFile <- modifiedExportMethylome(model=model$data, out.dir=out.dir, context=context, name=name)
    rm(model)
  }
  rm(methylome.data)
}

