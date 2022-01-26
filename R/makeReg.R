#' Generates cytosine clusters in genome
#'
#' Uses the extracted cytosines from CfromFASTA and finds clsuters of cytosines along the genome.
#'
#' @param ref.genome
#' @param contexts
#' @param chr
#' @param min.C
#' @param N.boot
#' @param N.sim.C
#' @param fp.rate
#' @param set.tol
#' @param out.dir
#' @param out.name
#' @param makeRegnull
#' @import GenomicRanges
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
#'
makeReg<-function(ref.genome, contexts, chr, min.C, N.boot, N.sim.C, fp.rate, set.tol, out.dir, out.name, makeRegnull)
{

  # Defining function for cleaning objects in the global environment from within the function
  CleanEnvir <- function(pattern)
  {
    rm(list = ls(envir=globalenv())[
      grep(pattern, ls(envir=globalenv()))], envir = globalenv())
  }

  #ref.genome<-out
  #context="CHH"
  #chr=1
  #min.C=5
  #N.boot=10^5
  #fp.rate=0.01
  #set.tol<-0.001

  # This function constructs the regions based on a sequence input
  conReg<-function(seq.in, min.C, win, chr)
  {

    start<-NULL
    stop<-NULL
    index<-rep(0, length(seq.in))

    pb <- txtProgressBar(min = 1, max = length(seq.in), char = "=", style = 3, file = "")
    for (i in seq_len(c(length(seq.in)- c(min.C-1))))
    {
      if (as.numeric(seq.in[i+ c(min.C-1)])-as.numeric(seq.in[i]) <= win){index[i:c(i + c(min.C-1))]<-rep(1, min.C)}

      Sys.sleep(1/length(seq.in)); setTxtProgressBar(pb, i)

    }

    close(pb)

    runs<-rle(index > 0)
    myruns = which(runs$values == TRUE & runs$lengths >= min.C)
    runs.lengths.cumsum = cumsum(runs$lengths)
    ends = runs.lengths.cumsum[myruns]
    newindex = ifelse(myruns>1, myruns-1, 0)
    starts = runs.lengths.cumsum[newindex] + 1
    if (0 %in% newindex) starts = c(1,starts)


    start.pos<-as.numeric(unlist(seq.in[starts]))
    end.pos<-as.numeric(unlist(seq.in[ends]))
    cluster.length<-end.pos - start.pos

    cluster<-cbind(chr, start.pos, end.pos, cluster.length)
    colnames(cluster)<-c("chr","start", "end", "cluster.length")
    cluster<-as.data.frame(cluster, stringsAsFactors = FALSE)
    cluster$chr<-as.character(cluster$chr)
    cluster$start<-as.integer(cluster$start)
    cluster$end<-as.integer(cluster$end)
    cluster$cluster.length<-as.integer(cluster$cluster.length)
    cluster$region<-paste("reg", 1:nrow(cluster), sep="")

    return(cluster)
  }




  cat("- Reading in the data", "\n")
  cat("-------------------------", "\n")
  cat("                         ", "\n")

  # Subselect genome based on chr.id
  #$ref.genome <- ref.genome %>% filter(ref.genome$chr==chr)

  # Making the position vectors
  CHH.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CHH"),]$pos
  CHG.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CHG"),]$pos
  CG.pos.plus<-ref.genome[which(ref.genome$strand  == "+" &  ref.genome$context == "CG"),]$pos
  CHH.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CHH"),]$pos
  CHG.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CHG"),]$pos
  CG.pos.minus<-ref.genome[which(ref.genome$strand  == "-" &  ref.genome$context == "CG"),]$pos

  # Determining the number of CHH, CHG and CG sites (+) strand
  N.CHH.obs<-length(CHH.pos.plus)
  N.CHG.obs<-length(CHG.pos.plus)
  N.CG.obs<-length(CG.pos.plus)


  if (N.sim.C == "all")
  {
    chr.length<-max(ref.genome$pos)
    N.CHH<-N.CHH.obs
    N.CHG<-N.CHG.obs
    N.CG<-N.CG.obs
    set.seed(123)
    rp.CHH<-sample(1:chr.length, size=N.CHH, replace = FALSE)
    rp.CHG<-sample(1:chr.length, size=N.CHG, replace = FALSE)
    rp.CG<-sample(1:chr.length, size=N.CG, replace = FALSE)
  }else{
    if (is.numeric(N.sim.C) == FALSE){stop("N.sim.C needs to be numeric")}else{
      chr.length<-max(ref.genome$pos)
      scale.factor<-N.sim.C/c(N.CHH.obs + N.CHG.obs + N.CG.obs)
      chr.length<-round(chr.length*scale.factor,0)
      N.CHH<-round(N.CHH.obs*scale.factor,0)
      N.CHG<-round(N.CHG.obs*scale.factor,0)
      N.CG<-round(N.CG.obs*scale.factor,0)
      set.seed(123)
      rp.CHH<-sample(1:chr.length, size=N.CHH, replace = FALSE)
      rp.CHG<-sample(1:chr.length, size=N.CHG, replace = FALSE)
      rp.CG<-sample(1:chr.length, size=N.CG, replace = FALSE)
    }
  }


  # Remove this large data. It's no longer needed
  CleanEnvir(pattern = "ref.genome")

  cat("- Simulating cytosines", "\n")
  cat("-------------------------", "\n")

  rp.CHH<-sort(rp.CHH, method="quick")
  rp.CHG<-sort(rp.CHG, method="quick")
  rp.CG<-sort(rp.CG, method="quick")

  CHH.lost<-0
  CHG.lost<-0
  CG.lost<-0

  tolerance = set.tol + 0.00001

  while (tolerance > set.tol)
  {
    ### CHH
    # Draw random positions (rp) of CHH sites
    rp.CHH<-c(rp.CHH, sample(1:chr.length, size= CHH.lost, replace =FALSE))
    if  (length(rp.CHH) > N.CHH){rp.CHH<-sample(rp.CHH, size = N.CHH, replace=F)}
    rp.CHH<-sort(rp.CHH)

    ### CHG
    # Draw random positions (rp) of CHG sites
    # Rules 1 to 3
    rp.CHG<-c(rp.CHG, sample(1:chr.length, size= CHG.lost, replace =FALSE))
    if  (length(rp.CHG) > N.CHG){rp.CHG<-sample(rp.CHG, size = N.CHG, replace=F)}
    rp.CHG<-sort(rp.CHG)
    rp.CHG<-rp.CHG[which(diff(rp.CHG, lag=1) >= 3)]

    # Rule 4
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, rp.CHG))]

    # Rule 5
    # no action

    # Rule 6
    # no action

    # Rule 7
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CHG+1)))]

    # Rule 8
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CHG+2)))]


    ### CG
    # Draw random positions (rp) of CG sites
    # Rules 1 to 3
    rp.CG<-c(rp.CG, sample(1:chr.length, size= CG.lost, replace =FALSE))
    if  (length(rp.CG) > N.CG){rp.CG<-sample(rp.CG, size = N.CG, replace=F)}
    rp.CG<-sort(rp.CG)
    rp.CG<-rp.CG[which(diff(rp.CG, lag=1) >= 2)]

    # Rule 4
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, rp.CG))]

    # Rule 5
    rp.CHH.temp<-rp.CHH+1
    rp.CHH.temp<-rp.CHH.temp[which(is.element(rp.CHH.temp, rp.CG))]
    rp.CHH<-rp.CHH[which(!is.element(c(rp.CHH+1), rp.CG))]
    rp.CHG<-c(rp.CHG, c(rp.CHH.temp-1))
    rp.CHG<-unique(rp.CHG)
    if (length(rp.CHG) > N.CHG){rp.CHG<-sample(rp.CHG, size = N.CHG, replace =F)}
    rp.CHG<-sort(rp.CHG)
    rp.CHG<-rp.CHG[which(diff(rp.CHG, lag=1) >= 3)]

    # Rule 6
    # no action

    # Rule 7
    rp.CHH<-rp.CHH[which(!is.element(rp.CHH, c(rp.CG+1)))]

    # Rule 8
    rp.CHG<-rp.CHG[which(!is.element(rp.CHG, rp.CG))]

    # Rule 9
    # no action

    # Rule 10
    rp.CHG.temp<-rp.CHG+2
    rp.CHG.temp<-rp.CHG.temp[which(is.element(rp.CHG.temp, rp.CG))]
    rp.CHG<-rp.CHG[which(!is.element(c(rp.CHG+2), rp.CG))]
    rp.CG<-c(rp.CG, c(rp.CHG.temp-2))
    rp.CG<-unique(rp.CG)
    if (length(rp.CG) > N.CG){rp.CG<-sample(rp.CG, size = N.CG, replace =F)}
    rp.CG<-sort(rp.CG)
    rp.CG<-rp.CG[which(diff(rp.CG, lag=1) >= 2)]

    # Rule 11
    rp.CHG<-rp.CHG[which(!is.element(rp.CHG, c(rp.CG+1)))]

    # Counting the lost sites
    CHH.lost<-N.CHH - length(rp.CHH)
    CHG.lost<-N.CHG - length(rp.CHG)
    CG.lost<-N.CG - length(rp.CG)
    N.lost<-CHH.lost + CHG.lost + CG.lost

    tolerance<-N.lost/(N.CG + N.CHG + N.CHH)

    cat(tolerance, "\n")

  }

  cat("                         ", "\n")
  cat("- Converged", "\n")
  cat("-------------------------", "\n")

  sim.out<-data.frame(c(N.CG.obs, N.CHG.obs, N.CHH.obs), c(length(rp.CG), length(rp.CHG), length(rp.CHH)))
  sim.out<-rbind(sim.out, colSums(sim.out))
  rownames(sim.out)<-c("CG", "CHG", "CHH", "Total")
  sim.out$percent<-sim.out[,2]/sim.out[,1]*100
  colnames(sim.out)<-c("N.observed", "N.simulated", "Percent")

  print(sim.out)

  cat("                         ", "\n")

  cat("- Calling regions:", "\n")
  cat("-------------------------", "\n")
  cat("                         ", "\n")


  for (icontext in 1:length(contexts))
  {
    context.temp<-contexts[icontext]


    if (context.temp == "C")
    {
      sim.geno.out<-c(rp.CG, rp.CHG, rp.CHH)
      sim.geno.out<-sort(sim.geno.out, method="quick")

      obs.geno.plus<-sort(c(CG.pos.plus, CHG.pos.plus, CHH.pos.plus))
      obs.geno.minus<-sort(c(CG.pos.minus, CHG.pos.minus, CHH.pos.minus))

      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] -sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))

      cat("Building from: C (+) strand", "\n")
      reg.obs.plus<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      cat("Building from: C (-) strand", "\n")
      reg.obs.minus<-conReg(seq.in=obs.geno.minus, min.C=min.C, win=null.dist, chr=chr)
      df <- data.frame(id=c(rep("+", nrow(reg.obs.plus)), rep("-", nrow(reg.obs.minus))),
                       start=c(reg.obs.plus[,2], reg.obs.minus[,2]), end=c(reg.obs.plus[,3], reg.obs.minus[,3]))
      gr <- GRanges(seqnames = rep(1,nrow(df)),IRanges(start = df$start, end = df$end))
      reg.obs<-as.data.frame(reduce(gr))[,1:3]
      reg.obs[,1]<-rep(chr, nrow(reg.obs))
      reg.obs$cluster.length <- reg.obs$end - reg.obs$start
      reg.obs$region<-paste("reg", 1:nrow(reg.obs), sep="")
      colnames(reg.obs)[1]<-"chr"

      if ("C" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: C null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
      }
      if (!is.element("C", contexts[which(makeRegnull == TRUE)]))
      {
        cat("C NULL omitted", "\n")
      }
    }

    if (context.temp == "CG")
    {
      sim.geno.out<-rp.CG

      obs.geno.plus<-sort(CG.pos.plus)
      obs.geno.minus<-sort(CG.pos.minus)

      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))

      cat("Building from: from CG (+) strand", "\n")
      reg.obs<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      reg.obs$start<-reg.obs$start - 1
      reg.obs$end<-reg.obs$end + 1
      reg.obs$cluster.length<-reg.obs$end - reg.obs$start

      if ("CG" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: CG null regions", "\n")
        cat("", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
        reg.sim$start<-reg.sim$start - 1
        reg.sim$end<-reg.sim$end + 1
        reg.sim$cluster.length<-reg.sim$end - reg.sim$start
      }
      if (!is.element("CG", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CG NULL omitted", "\n")
      }
    }

    if (context.temp == "CHG")
    {
      sim.geno.out<-rp.CHG
      obs.geno.plus<-sort(CHG.pos.plus)
      obs.geno.minus<-sort(CHG.pos.minus)

      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))

      cat("Building from: CHG (+) strand", "\n")
      reg.obs<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      reg.obs$start<-reg.obs$start - 2
      reg.obs$end<-reg.obs$end + 2
      reg.obs$cluster.length<-reg.obs$end - reg.obs$start

      if ("CHG" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("Building: CHG null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
        reg.sim$start<-reg.sim$start - 2
        reg.sim$end<-reg.sim$end + 2
        reg.sim$cluster.length<-reg.sim$end - reg.sim$start
      }
      if (!is.element("CHG", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CHG NULL omitted", "\n")
      }
    }

    if (context.temp == "CHH")
    {
      sim.geno.out<-rp.CHH
      obs.geno.plus<-sort(CHH.pos.plus)
      obs.geno.minus<-sort(CHH.pos.minus)

      set.seed(42)
      index1<-sample(1:c(length(sim.geno.out)-min.C), size=N.boot, replace=TRUE)
      index2<-index1 + c(min.C-1)
      dist<-sim.geno.out[index2] - sim.geno.out[index1]
      null.dist<-as.numeric(quantile(dist, probs=fp.rate))

      cat("Building from: CHH (+) strand", "\n")
      reg.obs.plus<-conReg(seq.in=obs.geno.plus, min.C=min.C, win=null.dist, chr=chr)
      cat("Building from: CHH (-) strand", "\n")
      reg.obs.minus<-conReg(seq.in=obs.geno.minus, min.C=min.C, win=null.dist, chr=chr)
      df <- data.frame(id=c(rep("+", nrow(reg.obs.plus)), rep("-", nrow(reg.obs.minus))),
                       start=c(reg.obs.plus[,2], reg.obs.minus[,2]), end=c(reg.obs.plus[,3], reg.obs.minus[,3]))
      gr <- GRanges(seqnames = rep(1,nrow(df)),IRanges(start = df$start, end = df$end))
      reg.obs<-as.data.frame(reduce(gr))[,1:3]
      reg.obs[,1]<-rep(chr, nrow(reg.obs))
      reg.obs$cluster.length <- reg.obs$end - reg.obs$start
      reg.obs$region<-paste("reg", 1:nrow(reg.obs), sep="")
      colnames(reg.obs)[1]<-"chr"


      if ("CHH" %in% contexts[which(makeRegnull == TRUE)])
      {
        cat("CHH null regions", "\n")
        reg.sim<-conReg(seq.in=sim.geno.out, min.C=min.C, win=null.dist, chr=chr)
      }
      if (!is.element("C", contexts[which(makeRegnull == TRUE)]))
      {
        cat("CHH NULL omitted", "\n")
      }
    }

    if (context.temp %in% contexts[which(makeRegnull == TRUE)])
    {
      output<-list(reg.obs, reg.sim, chr, context.temp)
      names(output)<-c("reg.obs", "reg.sim", "chr", "context")
      dput(output, paste(out.dir,"/",out.name,"_regions_","chr", chr, "_", context.temp, ".Rdata", sep=""))
    }
    if (!is.element(context.temp, contexts[which(makeRegnull == TRUE)]))
    {
      reg.sim<-NA
      output<-list(reg.obs, reg.sim, chr, context.temp)
      names(output)<-c("reg.obs", "reg.sim", "chr", "context")
      dput(output, paste(out.dir,"/",out.name, "_regions_","chr", chr, "_", context.temp, ".Rdata", sep=""))
    }

  } # End of context loop

  # Removing specific objects from global environemnt
  CleanEnvir(pattern = "reg.obs.plus")
  CleanEnvir(pattern = "reg.obs.minus")
  CleanEnvir(pattern = "obs.geno.plus")
  CleanEnvir(pattern = "obs.geno.minus")

}

