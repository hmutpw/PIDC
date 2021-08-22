#' Inferring Gene Regulatory Network (GRN) using partial information decomposition and context (PIDC) method
#'
#' The PIDC package implements the method in paper
#' (\href{https://linkinghub.elsevier.com/retrieve/pii/S2405471217303861}{2017 Chan et al.})
#' with few modifications. Besides, The code used for calculating the \code{Multivariate
#' Information} (MI) and \code{specific.information}, discretizing gene expression data
#' were adopted from package \code{\link{Informeasure}}.
#'
#' @param expMat A \code{matrix} storing the gene expression data. Rows corresponding
#' to features (eg. gene symbols) and columns corresponding to samples (cells).
#' Raw read counts, UMIs, TPMs or logNormalized counts were supported.
#' @param regulators The regulator genes used for GRN inferring (eg, transcription
#' factors). At least two regulators required. Default: NULL, using all the genes
#' in gene expression matrix.
#' @param targets The target genes used for GRN inferring. Default: NULL, using
#' all the genes in gene expression matrix.
#' @param logNormalized Whether the input have been logNormalized?
#' @param ncores Number of cores used for parallel calculation. The running time
#' heavily depend on the number of regulators and targets, for genes over 5000,
#' we strongly suggest you to use multi-cores. Default: 1.
#' @param diag Numeric. The weight in the diagonal of output square matrix. Default: 1.
#' @param verbose Whether to print message while running. Default: TRUE.
#'
#' @return A matrix with weighted value between regulators and targets.
#'
#' @importFrom purrr map map2 pmap transpose
#' @importFrom parallel makeCluster clusterExport clusterEvalQ stopCluster
#' @importFrom pbapply pbapply pblapply
#' @importFrom methods as is
#' @importFrom stats ecdf
#' @export
#'
#' @examples
#' data(expMat)
#' PIDC_res <- PIDC(expMat)
#' head(PIDC_res)
#'
#'
#'
PIDC <- function(expMat, regulators=NULL, targets=NULL, logNormalized=FALSE,
                 ncores=1, diag=c("auto","zero","one"), verbose = interactive()){
  if(verbose) message("[1] Filtering and Normalizing data...")
  .checkPIDCArgs(expMat=expMat, regulators=regulators, targets=targets,
                 logNormalized=logNormalized, ncores=ncores)
  diag <- match.arg(diag)
  #---check expMat
  expMat <- as.matrix(expMat)
  if(!logNormalized){
    expMat <- log2(expMat+1)
  }
  #---check regulators
  if(is.null(regulators)){
    regulators <- row.names(expMat)
  }
  regulators <- intersect(regulators, row.names(expMat))
  #---check targets
  if(is.null(targets)){
    targets <- row.names(expMat)
  }
  targets <- intersect(targets, row.names(expMat))
  ######calculating
  #1.discretize Gene
  if(verbose) message("[2] Discretizing expression matrix into bins...")
  expMat <- as.list(as.data.frame(t(expMat)))
  discret_list <- pbapply::pblapply(X=expMat,FUN = .discretizeGene)
  #2. Calculating the proportional unique contribution (PUC) score matrix
  if(verbose) message("[3] Calculating proportional unique contribution(PUC) ",
                      "matrix using ",length(regulators)," regulators and ",
                      length(targets)," targets...\n(This step may taken dozens ",
                      "of hours, please be patient)")
  #---parallel calculation
  cl <- parallel::makeCluster(ncores)
  parallel::clusterEvalQ(cl,library(purrr))
  parallel::clusterEvalQ(cl,library(dplyr))
  parallel::clusterEvalQ(cl,library(pbapply))
  PUC_list <- pbapply::pblapply(X = targets, FUN = .getPUC,
                                regulators = regulators,
                                discret_list = discret_list,
                                cl=cl)
  parallel::stopCluster(cl)
  names(PUC_list) <- targets
  PUC_mat <- do.call(cbind, PUC_list)
  #3. Summaring Uxy and calculating the weight matrix
  if(verbose) message("[4] Summaring Uxy and calculating the weight matrix...")
  out_res <- .FUxy(Uxy_mat = PUC_mat)
  if(diag=="auto"){
    return(out_res)
  }else if(diag=="zero"){
    diag(out_res) <- 0
  }else if(diag=="one"){
    diag(out_res) <- 1
  }
  out_res
}


#' Get regulatory network with PIDC matrix
#'
#' Inferring gene regulatory network using weighted square matrix from PIDC output.
#'
#' @param weightMat A square matrix whose element is positive values.
#' @param methods The name of the network inference algorithm. Default: aracne.
#' @param cutoff Set the cutoff of regulatory networks. Default: NULL
#'
#' @return A data.frame or a list with regulatory networks.
#' @importFrom minet aracne clr mrnet mrnetb
#' @importFrom reshape2 melt
#' @importFrom stats sd ks.test rnorm qnorm density dnorm
#' @importFrom utils installed.packages capture.output
#' @importFrom parallel makeCluster stopCluster
#' @importFrom pbapply pblapply
#'
#' @export
#'
#' @examples
#' data(expMat)
#' PIDC_res <- PIDC(expMat)
#' PIDC_net <- matToNet(PIDC_res)
matToNet <- function(weightMat,
                     methods = c("aracne", "clr", "mrnet", "mrnetb"),
                     cutoff = NULL,
                     ncores=1){
  methods <- match.arg(methods)
  mat <- as.matrix(weightMat)
  if(nrow(mat)!=ncol(mat)) stop("The input matrix must be square matrix!")
  if(min(mat)<0) stop("The input matrix can not have negative values!")
  if(!is.null(cutoff)){
    if(!is.numeric(cutoff) || !(cutoff>=0 && cutoff<=1) || length(cutoff)!=1){
      stop("The cutoff must be a numeric value between 0 and 1!")
    }
  }
  message("[1] Generating Gene regulatory networks using: ",methods,".")
  if(methods=="aracne"){
    net <- minet::aracne(mat)
  }else if(methods=="clr"){
    net <- minet::clr(mat)
  }else if(methods=="mrnet"){
    net <- minet::mrnet(mat)
  }else if(methods=="mrnetb"){
    net <- minet::mrnetb(mat)
  }
  out_tab <- reshape2::melt(data = net)
  colnames(out_tab) <- c("regulator","target","weight")
  if(is.null(cutoff)){
    message("[2] Inferring cutoff by estimateing distribution...")
    out_list <- split(x = out_tab, f = out_tab$regulator)
    cl <- parallel::makeCluster(ncores)
    out_flter_list <- pbapply::pblapply(X = out_list, FUN = .getWeightThreshold, cl=cl)
    parallel::stopCluster(cl)
    names(out_flter_list) <- names(out_list)
    out_res <- do.call(rbind,out_flter_list)
  }else{
    out_res <- out_tab[out_tab$weight>cutoff,]
  }
  row.names(out_res) <- NULL
  out_res[order(out_res$regulator,out_res$weight,out_res$target,
                decreasing = c(FALSE,TRUE,FALSE)),]
}

######
#---1. check input arguments
######
.checkPIDCArgs <- function(expMat, regulators, targets, logNormalized, ncores){
  #---check expMat
  if(!is.matrix(expMat) && !is.array(expMat) && !is(expMat,"dgCMatrix")){
    stop("The expMat must be a two-dimensional matrix where the row corresponds",
         " to a gene and each column corresponds to a sample")
  }
  if (length(dim(expMat)) != 2) {
    stop("The expMat must be a two-dimensional matrix where the row corresponds",
         " to a gene and each column corresponds to a sample")
  }
  if (is.null(rownames(expMat))) {
    stop("The expMat must contain the names of the genes as rownames.")
  }
  expMat <- as.matrix(expMat)
  if(!is.numeric(expMat)) stop("The expMat contain non-numeric values.")
  if(any(is.na(expMat))) stop("The expMat contain NA.")
  gene_names <- unique(rownames(expMat))
  #---check regulators in expMat
  regulators <- unique(regulators)
  if(!is.null(regulators)){
    if(!is.character(regulators)) stop("The regulators must characters!")
    if(length(intersect(gene_names, regulators))<2){
      stop("At least two regulators requiered, please check the names of the ",
           "expMat and the gene ids (names) in your regulators!")
    }else if(length(intersect(gene_names, regulators)) < length(regulators)){
      ratio <- length(setdiff(regulators,gene_names))/length(regulators)
      warning(length(setdiff(regulators,gene_names))," out of ",length(regulators),
              " (",round(ratio*100,2),"%) regulators not found in your expMat!")
    }
  }
  #---check targets in expMat
  targets <- unique(targets)
  if(!is.null(targets)){
    if(!is.character(targets)) stop("The targets must characters!")
    if(length(intersect(gene_names, targets))<1){
      stop("At least one targets requiered, please check the names of the ",
           "expMat and the gene ids (names) in your targets!")
    }else if(length(intersect(gene_names, targets)) < length(targets)){
      ratio <- length(setdiff(targets,gene_names))/length(targets)
      warning(length(setdiff(targets,gene_names))," out of ",length(targets),
              " (",round(ratio*100,2),"%) targets not found in your expMat!")
    }
  }
  if(!is.logical(logNormalized)) stop("The logNormalized must be TRUE(T) or FALSE(F)!")
  if(!is.numeric(ncores) || ncores<1) stop("The logNormalized must be positive numbers!")
}
######
#---2. discrete one gene expression into equal length items
######
.discretizeGene <- function(x, method = c("uniform_width","uniform_frequency")){
  method <- match.arg(method)
  if(method=="uniform_width"){
    numBins <- floor(sqrt(length(x)))
    r <- range(x)
    b <- seq(from = r[1], to = r[2], length.out = (numBins + 1))
    X <- cut(x, breaks = b , include.lowest = TRUE)
  }else if("uniform_frequency"){
    numBins <- floor(sqrt(length(x)))
    nrepl <- floor(length(x)/numBins)
    nplus <- sample(seq_len(numBins), length(x) - nrepl * numBins)
    nrep <- rep(nrepl, numBins)
    nrep[nplus] <- nrepl + 1
    X <- x[order(x)] <- rep(seq.int(numBins), nrep)
  }
  return(X)
}
######
#---3. get pairwise probability from  frequency table
######
# x discrete factor
# y discrete factor
.freqTable <- function(x, y, method=c("ML", "Jeffreys", "Laplace", "SG", "minimax")){
  method <- match.arg(method)
  if(is.list(x) || is.list(y)){
    x <- unlist(x)
    y <- unlist(y)
  }
  tab <- table(x,y)
  if(method=="ML"){
    probs <- tab/sum(tab)
  }else if(method=="Jeffreys"){
    probs <- (tab + 0.5)/sum(tab + 0.5)
  }else if(method=="Laplace"){
    probs <- (tab + 1)/sum(tab + 1)
  }else if(method=="SG"){
    probs <- (tab + length(tab))/sum(tab + length(tab))
  }else if(method=="minimax"){
    probs <- (tab + sqrt(sum(tab))/length(tab))/sum(tab + sqrt(sum(tab))/length(tab))
  }
  probs
}

######
#---3.get all regulator and one target pair PUC scores
######
.getPUC <- function(target, discret_list, regulators = names(discret_list)){
  if(length(regulators)<2) stop("At least two regulators needed!")
  if(length(target)!=1) stop("Only support one target each time!")
  if(is.null(names(discret_list))) stop("The names of discret_list must be gene symbols!")
  Z <- discret_list[target]
  Xi <- discret_list[regulators]
  p_XiZ <- purrr::map2(.x = Xi, .y = Z, .f = .freqTable)
  p_Xi <- purrr::map(.x = p_XiZ, .f = rowSums)
  p_Z <- purrr::map(.x = p_XiZ, .f = colSums)
  I_XiZ <- purrr::pmap(.l=list(freqs=p_XiZ, p_i=p_Xi, p_z=p_Z), .f = .MI)
  Ispec_XiZ <- purrr::pmap(.l=list(p_iz=p_XiZ, p_i=p_Xi, p_z=p_Z), .f = .specific.information)
  U_xy <- .puc_per_target(I_XiZ = I_XiZ, Ispec_XiZ = Ispec_XiZ, p_Z = p_Z)
  return(U_xy)
}

#' I_XiZ, Ispec_XiZ and p_Z is list with equally length, theoretically, p_Z is the same value.
.puc_per_target <- function(I_XiZ, Ispec_XiZ, p_Z){
  gene_names <- names(Ispec_XiZ)
  Ispec_XiZ_bins <- lapply(purrr::transpose(Ispec_XiZ),unlist)
  Ispec_XiZ_bins_pos <- which(sapply(Ispec_XiZ_bins,sum)>0)
  p_Z_bins <- lapply(purrr::transpose(p_Z),unlist)
  p_Z_bins_pos <- which(sapply(p_Z_bins,sum)>0)
  if(!identical(Ispec_XiZ_bins_pos, p_Z_bins_pos)){
    warning("The bins with values may not match between Ispec and p_Z!")
  }
  overlap_bins_pos <- union(Ispec_XiZ_bins_pos, p_Z_bins_pos)
  Ispec_XiZ_bins_filter <- Ispec_XiZ_bins[overlap_bins_pos]
  p_Z_bins_filter <- p_Z_bins[overlap_bins_pos]
  #---calculate sum redundance of Xi Z for all Yi
  per_bin_redundancy <- purrr::map2(.x = Ispec_XiZ_bins_filter,
                                    .y = p_Z_bins_filter,
                                    .f = function(x, y){
                                      if(is.null(names(x)) || is.null(names(y))){
                                        stop("The Ispec and p_Z in each bin must have the same gene names!")
                                      }
                                      x_sort <- sort(x)
                                      n_gene <- length(x_sort)
                                      x_cumsum <- cumsum(x_sort)
                                      x_over_num <- n_gene-order(x_sort)-1
                                      names(x_over_num) <- names(x_sort)
                                      x_over_sum <- x_sort*x_over_num
                                      x_value <- x_cumsum+x_over_sum
                                      y_value <- y[names(x_value)]
                                      out <- x_value*y_value
                                      out[names(x)]
                                    })
  #---check gene name
  check_gnames <- sapply(per_bin_redundancy,function(x,gname){
    identical(names(x),gname)},gname=gene_names)
  if(!all(check_gnames)){
    per_bin_redundancy <- lapply(per_bin_redundancy,function(x,gnames){
      x[gnames]},gname=gene_names)
  }
  per_gene_redundancy <- sapply(purrr::transpose(per_bin_redundancy),
                                function(x){sum(unlist(x))})
  per_gene_puc <- length(gene_names)-1-per_gene_redundancy/unlist(I_XiZ)[names(per_gene_redundancy)]
  return(per_gene_puc)
}

#---calculation of specific.information for each regulator-target pair
.specific.information <- function(p_iz, p_i, p_z, unit = c("log", "log2", "log10")){
  unit <- match.arg(unit)
  p_i_z <- t(t(p_iz)/p_z)  ##p(i|z)
  p_z_i <- t(p_iz/p_i) ##p(z|i)
  tmp <- t((log(1/p_z)) - (log(1/p_z_i))) * p_i_z
  if (unit == "log2")  tmp <- t((log(1/p_z, 2))  - (log(1/p_z_i, 2))) * p_i_z  # change from log to log2 scale
  if (unit == "log10") tmp <- t((log(1/p_z, 10)) - (log(1/p_z_i, 10))) * p_i_z # change from log to log10 scale
  colSums(tmp, na.rm = TRUE)
}

#---calculation of mutual information
.MI <- function(freqs, p_i, p_z, unit = c("log", "log2", "log10")){
  if(length(dimnames(freqs))!=2) stop("The dims of freqs must be 2!")
  unit <- match.arg(unit)
  MI <- .H(p_i, unit = unit) + .H(p_z, unit = unit) - .H(freqs, unit = unit)
  MI
}

.H <- function(freqs, unit = c("log", "log2", "log10")){
  unit = match.arg(unit)
  #freqs = freqs/sum(freqs)
  H = -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
  if (unit == "log2") H = H/log(2)
  if (unit == "log10") H = H/log(10)
  H
}

######
#---4. get the cumulative distribution score of PUC
######
.FUxy <- function(Uxy_mat){
  mat <- as.matrix(Uxy_mat)
  if(nrow(mat)==ncol(mat)){
    mat <- mat + t(mat)
  }
  F_x <- t(apply(mat,1,function(x){stats::ecdf(x)(x)}))
  colnames(F_x) <- colnames(mat)
  F_y <- apply(mat,2,function(x){stats::ecdf(x)(x)})
  row.names(F_y) <- row.names(mat)
  F_xy <- (F_x+F_y)*0.5
  return(F_xy)
}

######
#---get Weight Threshold
######

#---calculate the cutoff of score using method from AUCell package with some modification
.getWeightThreshold <- function(df, smallestPopPercent=.1,
                                densAdjust=2, thrP=0.01, nBreaks=100){
  #---progress bar
  if(length(unique(df[["regulator"]]))!=1) stop("Only support one regulator input!")
  auc <- df[["weight"]]
  names(auc) <- df[["target"]]
  gSetName <- unique(df[["regulator"]])

  nCells <- length(auc)
  skipGlobal <- TRUE
  skipRed <- FALSE
  skipSmallDens <- FALSE
  commentMsg <- ""
  aucThrs <- c()

  notPopPercent <- 1 - smallestPopPercent
  if(sum(auc==0) > (nCells*notPopPercent))
  {
    skipGlobal <- FALSE
    commentMsg <- paste(commentMsg,
                        round((sum(auc==0)/nCells)*100),
                        "% (more than ", notPopPercent,"%) of AUC are zero. ", sep="")
  }

  meanAUC <- mean(auc)
  sdAUC <- sd(auc)
  maybeNormalDistr <- !suppressWarnings(
    ks.test(auc, rnorm(max(100,length(auc)),mean=meanAUC, sd = sdAUC),
            alternative = "less")$p.value < .01)
  if(maybeNormalDistr){
    commentMsg <- paste0(commentMsg,
                         "The AUC might follow a normal distribution (random gene-set?). ")
    skipGlobal <- FALSE

    # aucThrs["outlierOfGlobal"] <- meanAUC + 2*sdAUC
    aucThrs["outlierOfGlobal"] <- qnorm(1-(thrP/nCells), mean=meanAUC, sd=sdAUC)
  }

  #V6
  histogram <- hist(c(0, auc/max(auc)), breaks=100, plot=FALSE)$count
  if((sum(histogram[1:5]) / sum(histogram)) >= notPopPercent*.75) {
    skipGlobal <- FALSE
    skipRed <- TRUE
    skipSmallDens <- TRUE
  }
  if((sum(histogram[1:10]) / sum(histogram)) >= notPopPercent*.50) {
    skipSmallDens <- TRUE
    skipGlobal <- FALSE
    # skipRed <- TRUE ?
    aucThrs["tenPercentOfMax"] <- max(auc)*.10
  }
  # print(skipRed)

  densCurve <- density(auc, adjust=densAdjust, cut=0)
  maximumsDens <- NULL
  inflPoints <- diff(sign(diff(densCurve$y)))
  maximumsDens <- which(inflPoints==-2)
  globalMax <- maximumsDens[which.max(densCurve$y[maximumsDens])]
  minimumDens <- which(inflPoints==2)
  smallMin <- NULL
  if(!skipSmallDens)
    smallMin <- data.table::last(minimumDens[which(minimumDens < globalMax)]) #1prev to max
  minimumDens <- c(smallMin,
                   minimumDens[which(minimumDens > globalMax)]) # all after maximum

  # Density-based threshold (V4):
  # First minimum after the biggest maximum   (adjust=2)
  densTrh <- NULL
  if(length(minimumDens)>0) # && (!skipMinimumDens))
  {
    densTrh <- densCurve$x[min(minimumDens)]
    # Commented on V6
    # Only keep if it is a real inflextion point
    # (i.e. next max at least 5% of the global max)
    if(length(maximumsDens)>0)
    {
      nextMaxs <- maximumsDens[which(densCurve$x[maximumsDens] > densTrh)]
      if((max(densCurve$y[nextMaxs])/max(densCurve$y))<.05)
      {
        densTrh <- NULL
        # print(gSetName)
      }
    }
  }

  ## TO DO: Check special cases with many zeroes
  auc <- sort(auc)
  distrs <- list()
  distrs[["Global_k1"]] <- list(mu=c(meanAUC, NA), sigma=c(sdAUC, NA), x=auc)


  if("mixtools" %in% rownames(installed.packages()))
  {
    na <- capture.output(distrs[["k2"]] <-
                           tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=2, verb=FALSE),
                                    # With fast, if there are many zeroes, it fails quite often
                                    error = function(e) {
                                      return(NULL)
                                    }))

    na <- capture.output(distrs[["k3"]] <-
                           tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=3, verb=FALSE),
                                    error = function(e) {
                                      return(NULL)
                                    }))

    if(is.null(distrs[["k2"]]) && is.null(distrs[["k3"]]))
    {
      if(sum(auc==0)<(nCells*notPopPercent*.5))
        skipGlobal <- FALSE    # only if not too many zeroes??
    }

    if(!is.null(distrs[["k2"]]))
    {
      compL <- which.min(distrs[["k2"]][["mu"]])
      compR <- which.max(distrs[["k2"]][["mu"]])
      ### Check distributions
      # Second distribution is "taller" than first one
      height1 <- .4/distrs[["k2"]][["sigma"]][compL]*
        distrs[["k2"]][["lambda"]][compL]
      height2 <- .4/distrs[["k2"]][["sigma"]][compR]*
        distrs[["k2"]][["lambda"]][compR]
      taller <- height1 < height2
      # Use global distr:
      # Mean of the global distr is included within the SD of the first
      # & Both means are included within the mean+SD of the Global distribution
      globalInclInFirst <-
        (distrs[["Global_k1"]]$mu[1] <
           (distrs[["k2"]][["mu"]][compL]+(1.5*distrs[["k2"]][["sigma"]][compL])))
      includedInGlobal <-
        ((distrs[["k2"]][["mu"]][compL] >
            (distrs[["Global_k1"]]$mu[1]-distrs[["Global_k1"]]$sigma[1])) &&
           (distrs[["k2"]][["mu"]][compR] <
              (distrs[["Global_k1"]]$mu[1]+distrs[["Global_k1"]]$sigma[1])))
      if(taller || (globalInclInFirst && includedInGlobal))
      {
        skipGlobal <- FALSE

        if(globalInclInFirst && includedInGlobal)
          commentMsg <- paste(commentMsg,
                              "The global distribution overlaps the partial distributions. ")
        if(taller && !includedInGlobal)
          commentMsg <- paste(commentMsg, "The right distribution is taller. ")
      }
    }
  }else{
    warning("Package 'mixtools' is not available to calculate the sub-distributions.")
  }

  glProb <- 1-(thrP/nCells + smallestPopPercent)   ## CORRECT?!?!
  aucThrs["Global_k1"] <- qnorm(glProb,# qnorm(1-(thrP/nCells),
                                mean=distrs[["Global_k1"]][["mu"]][1],
                                sd=distrs[["Global_k1"]][["sigma"]][1])
  if(!is.null(distrs[["k2"]]))
  {
    k2_L <- which.min(distrs[["k2"]][["mu"]]) # (sometimes the indexes are shifted)
    aucThrs["L_k2"] <- qnorm(1-(thrP/nCells),
                             mean=distrs[["k2"]][["mu"]][k2_L],
                             sd=distrs[["k2"]][["sigma"]][k2_L])
  }

  if(!is.null(distrs[["k3"]]))
  {
    k3_R <- which.max(distrs[["k3"]][["mu"]]) # R: right distribution
    k3_R_threshold <- qnorm(thrP,
                            mean=distrs[["k3"]][["mu"]][k3_R],
                            sd=distrs[["k3"]][["sigma"]][k3_R])
    if(k3_R_threshold > 0) aucThrs["R_k3"] <- k3_R_threshold
  }

  if(!is.null(densTrh))
  {
    aucThrs["minimumDens"] <- densTrh
  }

  aucThr <- aucThrs
  if(skipGlobal)
    aucThr <- aucThrs[which(!names(aucThrs) %in% "Global_k1")]
  # TO DO: Decide when to merge with GLOBAL

  if(skipRed)
    aucThr <- aucThrs[which(!names(aucThrs) %in% "L_k2")]
  # TO DO: Decide when to merge with GLOBAL

  aucThr <- aucThr[which.max(aucThr)] # to keep name
  if((length(aucThr)>0) && (names(aucThr) == "minimumDens"))
  {
    maximumsDens <- maximumsDens[which(densCurve$y[maximumsDens]>1)]
    if(length(maximumsDens) > 2)
    {
      tmp <- cbind(minimumDens[seq_len(length(maximumsDens)-1)],
                   maximumsDens[-1])
      FCs <- densCurve$y[tmp[,2]]/densCurve$y[tmp[,1]]
      if(any(FCs > 1.5))
        warning(gSetName,
                ":\tCheck the AUC histogram. ",
                "'minimumDens' was selected as the best threshold, ",
                "but there might be several distributions in the AUC.")
    }
  }

  if("minimumDens" %in% names(aucThrs))
    aucThr <- aucThrs["minimumDens"]
  if(length(aucThr)==0) aucThr <-  aucThrs[which.max(aucThrs)]
  df[df$weight>=aucThr,]
}




