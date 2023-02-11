#https://github.com/mrhelmus/phylogeny_manipulation/blob/master/AIC_func.r
AICc.phylolm<-function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
  
  if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
  LL <- logLik(mod)$logLik
  K <- logLik(mod)$df  #extract correct number of parameters included in model - this includes LM
  if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
  if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
  return(AICc)
}

# #https://stackoverflow.com/questions/44737985/r-remove-outliers-in-a-dataframe-grouped-by-factor
# remove_outliers <- function(x, na.rm = TRUE, ...) {
#   qnt <- quantile(x, probs=c(.1, .9), na.rm = na.rm, ...)
#   H <- 1.5 * IQR(x, na.rm = na.rm)
#   y <- x
#   y[x < (qnt[1] - H)] <- NA
#   y[x > (qnt[2] + H)] <- NA
#   y
# }

#https://stackoverflow.com/questions/44737985/r-remove-outliers-in-a-dataframe-grouped-by-factor
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.2, 0.8), na.rm = na.rm, ...)
  #H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1])] <- NA
  y[x > (qnt[2])] <- NA
  y
}


#function to match nodes from consensus 
#to individual gene trees with uneven sampling
#derived from Liam Revell's example-- need to test
match_nodes<-function(t1, t2){
  ## step one drop tips
  t1p<-drop.tip(t1,setdiff(t1$tip.label, t2$tip.label))
  t2p<-drop.tip(t2,setdiff(t2 $tip.label, t1$tip.label))
  
  ## step two match nodes "descendants"
  M<-matchNodes(t1p,t2p)
  
  ## step two match nodes "distances"
  M1<-matchNodes(t1,t1p,"distances")
  M2<-matchNodes(t2,t2p,"distances")
  
  ## final step, reconcile
  MM<-matrix(NA,t1$Nnode,2,dimnames=list(NULL,c("left","right")))
  
  for(i in 1:nrow(MM)){
    MM[i,1]<-M1[i,1]
    nn<-M[which(M[,1]==M1[i,2]),2]
    if(length(nn)>0){	
      if(length(which(M2[,2]==nn))>0){
        MM[i,2]<-M2[which(M2[,2]==nn),1]
      }
    } else {
    }   
  }
  return(MM)	
}


#function for matching edges between a phylogram and a cladeogram and estimating rates
#http://blog.phytools.org/2017/12/matching-edges-between-topologically.html
#trees is a list of trees where the first tree is a cladogram and the second tree
#is a phylogram, just by convention (does not check)
get_matched_edgelengths <- function(trees){
  M1<-matrix(NA,Ntip(trees[[1]]),length(trees),
             dimnames=list(trees[[1]]$tip.label,
                           paste("t[[",1:length(trees),"]]",sep="")))
  M2<-matrix(NA,trees[[1]]$Nnode,length(trees),
             dimnames=list(1:trees[[1]]$Nnode+Ntip(trees[[1]]),
                           paste("t[[",1:length(trees),"]]",sep="")))
  for(i in 1:length(trees)){
    M1[,i]<-phytools::matchLabels(trees[[1]],trees[[i]])[,2]
    M2[,i]<-phytools::matchNodes(trees[[1]],trees[[i]])[,2]
  }
  M<-rbind(M1,M2)
  #print(M)
  M<-M[-Ntip(trees[[1]])-1,] ## trim root node from M
  E<-matrix(NA,nrow(M),ncol(M),dimnames=dimnames(M))
  for(i in 1:ncol(M)) for(j in 1:nrow(M))
    E[j,i]<-trees[[i]]$edge.length[which(trees[[i]]$edge[,2]==M[j,i])]
  
  #names<-rownames(E)
  #distance/time = rate
  #E <- setNames(E[,2]/E[,1],names)
  
  return(E)
  
}

#https://rdrr.io/github/bozenne/butils/src/R/model.frame.gls.R
model.matrix.gls <- function(object, ...)
  model.matrix(terms(object), data = getData(object), ...)

model.frame.gls <- function(object, ...)
  model.frame(formula(object), data = getData(object), ...)

terms.gls <- function(object, ...)
  terms(model.frame(object),...)

#https://rdrr.io/github/bozenne/butils/src/R/model.frame.gls.R
model.matrix.phylolm <- function(object, ...)
  model.matrix(terms(object), data = getData(object), ...)

model.frame.phylolm <- function(object, ...)
  model.frame(formula(object), data = getData(object), ...)

terms.phylolm <- function(object, ...)
  terms(model.frame(object),...)



#https://github.com/rsquaredacademy/olsrr/blob/master/R/ols-added-variable-plot.R
{
ols_plot_added_variable <- function(model, print_plot = TRUE) {
  
  #check_model(model)
  
  data    <- ols_prep_avplot_data(model)
  xnames  <- colnames(data)
  nl      <- length(xnames)
  myplots <- list()
  
  for (i in 2:nl) {
    
    x <- ols_prep_regress_x(data, i)
    y <- ols_prep_regress_y(data, i)
    d <- data.frame(x, y)
    
    p <-
      eval(
        substitute(
          ggplot(d, aes(x = x, y = y)) +
            geom_point(colour = "blue", size = 2) +
            stat_smooth(method = "lm", se = FALSE) +
            xlab(paste(xnames[i], " | Others")) +
            ylab(paste(xnames[1], " | Others")),
          list(i = i)
        )
      )
    
    j <- i - 1
    myplots[[j]] <- p
    
  }
  
  if (print_plot) {
    marrangeGrob(myplots, nrow = 2, ncol = 2, top = "Added Variable Plots")
  } else {
    return(myplots)
  }
  
}

ols_prep_regress_x <- function(data, i) {
  
  x <- remove_columns(data, i)
  y <- select_columns(data, i)
  #lsfit(x, y)$residuals
  gls(y ~ x, correlation = corBrownian(phy = prunetree.shallow.time), method = "ML")$residuals
  
}

ols_prep_regress_y <- function(data, i) {
  
  x <- remove_columns(data, i)
  y <- select_columns(data)
  #lsfit(x, y)$residuals
  gls(y ~ x, correlation = corBrownian(phy = prunetree.shallow.time), method = "ML")$residuals
  
}

remove_columns <- function(data, i) {
  as.matrix(data[, c(-1, -i)])
}

select_columns <- function(data, i = 1) {
  as.matrix(data[, i])
}
}


gls.ci<-function (Y, X, Sigma) 
{
  n <- length(X)
  tr <- sum(diag(Sigma))
  Sigma <- n * Sigma/tr
  invSigma <- solve(Sigma)
  X1 <- rep(1, n)
  q <- 2
  C1 <- solve(t(X1) %*% invSigma %*% X1)
  Y_PGLSmean <- c(C1 %*% t(X1) %*% invSigma %*% Y)
  Y_PGLSdeviations = Y - Y_PGLSmean
  Y_PGLSvariance = (t(Y_PGLSdeviations) %*% invSigma %*% Y_PGLSdeviations)/(n - 
                                                                              1)
  SE_Y_mean = sqrt(Y_PGLSvariance/n)
  XX <- cbind(rep(1, n), X)
  C <- solve(t(XX) %*% invSigma %*% XX)
  w <- C %*% t(XX) %*% invSigma
  B <- w %*% Y
  Yhat <- XX %*% B
  Yresid = Y - Yhat
  Y_MSEresid <- c((t(Yresid) %*% invSigma %*% Yresid)/(n - 
                                                         q))
  a <- B[1]
  b <- B[2]
  SEa <- sqrt(diag(C) * Y_MSEresid)[1]
  SEb <- sqrt(diag(C) * Y_MSEresid)[2]
  intercept <- cbind(a, SEa)
  slope <- cbind(b, SEb)
  model <- rbind(intercept, slope)
  colnames(model) <- c("Estimate", "Std.Error")
  rownames(model) <- c("intercept", "slope")
  SEYhat <- sqrt(diag(XX %*% C %*% t(XX)) %*% ((t(Yresid) %*% 
                                                  invSigma %*% Yresid)/(n - q)))
  CI <- cbind(X, Yhat, SEYhat)
  Lower2.5 <- Yhat - qt(0.9, n) * SEYhat
  Lower5 <- Yhat - qt(0.95, n) * SEYhat
  Upper5 <- Yhat + qt(0.95, n) * SEYhat
  Upper2.5 <- Yhat + qt(0.9, n) * SEYhat
  CI <- cbind(CI, Lower2.5)
  CI <- cbind(CI, Lower5)
  CI <- cbind(CI, Upper5)
  CI <- cbind(CI, Upper2.5)
  CI <- CI[order(CI[, 1]), ]
  colnames(CI) <- c("X", "Yhat", "SEYhat", "Lower2.5", "Lower5", 
                    "Upper5", "Upper2.5")
  CI <- as.data.frame(CI)
  Xi <- seq(c(min(X) - abs(max(X))), to = c(abs(max(X)) * 5), 
            length.out = 100)
  Z <- c(Xi)
  ZZ <- cbind(rep(1, length(Z)), Z)
  SEYhat.Xi <- sqrt(diag(ZZ %*% C %*% t(ZZ)) %*% ((t(Yresid) %*% 
                                                     invSigma %*% Yresid)/(n - q)))
  Yhat.Xi <- a + b * Xi
  Lower2.5.Yhat.Xi <- Yhat.Xi - qt(0.9, n) * SEYhat.Xi
  Lower5.Yhat.Xi <- Yhat.Xi - qt(0.95, n) * SEYhat.Xi
  Upper5.Yhat.Xi <- Yhat.Xi + qt(0.95, n) * SEYhat.Xi
  Upper2.5.Yhat.Xi <- Yhat.Xi + qt(0.9, n) * SEYhat.Xi
  CI.plot <- cbind(Xi, Yhat.Xi, SEYhat.Xi, Lower2.5.Yhat.Xi, 
                   Lower5.Yhat.Xi, Upper5.Yhat.Xi, Upper2.5.Yhat.Xi)
  colnames(CI.plot) <- c("X", "Yhat", "SEYhat", "Lower2.5", 
                         "Lower5", "Upper5", "Upper2.5")
  CI.plot <- as.data.frame(CI.plot)
  results <- list(model, CI, CI.plot)
  names(results) <- c("model", "CI", "CI.plot")
  return(results)
}




#identify descendant node number for a given edge number
node_indices_edge<-function(tree, edges){
  edgetable<-tree$edge
  result<-list()
  
  for(i in 1:length(edges)){
    result[i]<- edgetable[edges[i],][2]
  }
  return(unlist(result))
}



anova.pgls.fixed <- function (object) 
{
  data <- object$data
  tlabels <- attr(terms(object$formula), "term.labels")
  k <- object$k
  n <- object$n
  NR <- length(tlabels) + 1
  rss <- resdf <- rep(NA, NR)
  rss[1] <- object$NSSQ
  resdf[1] <- n - 1
  lm <- object$param["lambda"]
  dl <- object$param["delta"]
  kp <- object$param["kappa"]
  for (i in 1:length(tlabels)) {
    fmla <- as.formula(paste(object$namey, " ~ ", paste(tlabels[1:i], collapse = "+")))
    plm <- pgls(fmla, data, lambda = lm, delta = dl, kappa = kp)
    rss[i + 1] <- plm$RSSQ
    resdf[i + 1] <- (n - 1) - plm$k + 1
  }
  ss <- c(abs(diff(rss)), object$RSSQ)
  df <- c(abs(diff(resdf)), n - k)
  ms <- ss/df
  fval <- ms/ms[NR]
  P <- pf(fval, df, df[NR], lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f = fval, P)
  table[length(P), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Residuals"), c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  structure(table, heading = c("Analysis of Variance Table", sprintf("Sequential SS for pgls: lambda = %0.2f, delta = %0.2f, kappa = %0.2f\n", lm, dl, kp), paste("Response:", deparse(formula(object)[[2L]]))), class = c("anova", "data.frame"))
}


#functions for calculating DR
DR_statistic <- function(x, return.mean = FALSE){
  
  rootnode <- length(x$tip.label) + 1
  
  sprates <- numeric(length(x$tip.label))
  for (i in 1:length(sprates)){
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode){
      el <- x$edge.length[x$edge[,2] == node]
      node <- x$edge[,1][x$edge[,2] == node]
      
      qx <- qx + el* (1 / 2^(index-1))
      
      index <- index + 1
    }
    sprates[i] <- 1/qx
  }
  
  if (return.mean){
    return(mean(sprates))		
  }else{
    names(sprates) <- x$tip.label
    return(sprates)
  }
  
}


z_transform<-function(data){
(data - mean(data))/sd(data)
}

# Interface to traitRate program (Mayrose & Otto, 2011)
# PACKAGE: ips
# CALLED BY: user
# AUTHOR: Christoph Heibl
# LAST UPDATE: 2014-08-07

traitRate <- function(phy, seq, x, mainType = "Optimize_Model", n,
                      charModelParam1 = 0.5, charModelParam2 = 1, 
                      gammaParam = 0.5, seqModelParam1 = 2,
                      exec = "/Applications/traitRate-1.1/programs/traitRate"){
  
  ## check input data
  ## ----------------
  if ( !inherits(phy, "phylo") )
    stop("'phy' is not of class 'phylo'")
  if ( !is.ultrametric(phy) )
    stop("'phy' must be ultrametric")
  phy$node.label <- NULL # traitRate does not parse node labels
  if ( !inherits(seq, "DNAbin") )
    stop("'seq' is not of class 'DNAbin'")
  if ( !is.matrix(seq) )
    stop("'seq' must be a matrix")
  if ( ncol(seq) > 15000 )
    stop("traitRate cannot handle > 15000 bp")
  
  ## write input data
  ## ----------------
  fn <- c("model.tree", "in_msa.fasta", "in_chars.fasta")
  write.tree(phy, fn[1])
  write.fas(seq, fn[2])
  write.fas(x, fn[3])
  
  outDir <- "RESULTS"
  outFile <- "traitRate.res"
  
  ## what type of analysis?
  ## ----------------------
  mainType <- match.arg(mainType, 
                        c("Optimize_Model", 
                          "runTraitBootstrap"))
  
  ## parametric bootstrapping
  ## ------------------------
  if ( mainType == "runTraitBootstrap" ){
    res.fn <- paste(outDir, outFile, sep = "/")
    if ( !file.exists(res.fn) ) stop("cannot find rate estimates")
    ## parse estimates of rates of trait evolution
    r.est <- traitRateOutput(res.fn, "rates")
    charModelParam1 <- r.est[1]
    charModelParam2 <- r.est[2]
    ll <- traitRateOutput(res.fn, "likelihoods")
    outDir <- "RESULTS_BOOTSTRAP"
  } else {
    out <- NULL
  }
  ## write parameters file
  ## ---------------------
  traitRateParams(mainType = mainType, n, fn, 
                  outDir, outFile,
                  charModelParam1 = charModelParam1, 
                  charModelParam2 = charModelParam2,
                  gammaParam, seqModelParam1)
  
  call <- paste(exec, "traitRate.doubleRep", sep = "/")
  call <- paste(call, "params.txt")
  system(call)
  
  ## run traitRate on simulated replicates
  ## -------------------------------------
  if ( mainType == "runTraitBootstrap" ){
    for ( i in seq(from = 0, to = n - 1) ){
      od <- paste(outDir, "/sim_", i, sep = "")
      fn[3] <- paste(od, "simRandomChars.fasta", sep = "/")
      traitRateParams(mainType = "Optimize_Model", n, fn, 
                      outDir = od, outFile = outFile,
                      charModelParam1, charModelParam2,
                      gammaParam, seqModelParam1)
      
      call <- paste(exec, "traitRate.doubleRep", sep = "/")
      call <- paste(call, "params.txt")
      system(call)
      
    }
    res.fn <- paste("sim", 0:(n - 1), sep = "_")
    res.fn <- paste(outDir, res.fn, outFile, 
                    sep = "/")
    out <- lapply(res.fn, traitRateOutput)
    out <- do.call(rbind, out)
    out <- cbind(out, 
                 diff = out[, "logL"] - out[, "logL0"])
  }
  out
}


traitRateParams <- function(mainType, n, fn, outDir, outFile,
                            charModelParam1, charModelParam2,
                            gammaParam, seqModelParam1){
  
  ## assemble control file
  ## ---------------------
  ctrl <- c(paste("_mainType", mainType),
            paste("_treeFile", fn[1]),
            paste("_characterFile", fn[3]),
            paste("_seqFile", fn[2]),
            paste("_outDir", outDir),
            paste("_outFile", outFile),
            "_logFile log.txt",
            "_scaledTreeFile scaled.tree",
            paste("_charModelParam1", charModelParam1),
            paste("_charModelParam2", charModelParam2),
            paste("_gammaParam", gammaParam),
            paste("_seqModelParam1", seqModelParam1),
            "_relRate 1",
            "_seqModelType HKY",
            "_logValue 3",
            "_bScaleTree 1",
            "_stochasicMappingIterations 100",
            "_treeLength 1.0")
  
  ## number of iterations
  if ( mainType == "runTraitBootstrap" ){
    if ( missing(n) ) n <- 200
    ctrl <- c(ctrl,
              paste("_", c("start", "end"),
                    "SimulationsIter ", c(0, n - 1), sep = ""))
  }
  write(ctrl, file = "params.txt")
}


traitRateOutput <- function(file, what = "likelihoods"){
  
  what <- match.arg(what, c("rates", "likelihoods"))
  sep <- ifelse(what == "rates", "", "\n")
  res <- scan(file, what = "c", sep = sep, 
              quiet = TRUE)
  if ( what == "rates" ){
    res <- res[grep("charModelParam", res)]
    res <- as.numeric(gsub("^.+=", "", res))
    names(res) <- c("r01", "r10")
  }
  if ( what == "likelihoods" ){
    id <- grep("LogLikelihood Model 0", res)
    res <- as.numeric(gsub("^.+= ", "", res[-1:0 + id]))
    names(res) <- c("logL", "logL0")
  }
  res
}


#https://raw.githubusercontent.com/evolucionario/fossilgraft/master/fossil.graft
fossil.graft <- function(phy, tip, fossil, fossil.age, edge.rel.length=0.5) {
  
  # Identify the edge where the new branch will be attached to
  
  if(length(tip) == 1) {
    tip0 <- which(phy$tip.label == tip) 
    mrcaage0 <- 0
  }
  else {
    tip0 <- getMRCA(phy, tip=tip)
    mrcaage0 <- branching.times(phy)[as.character(tip0)]
  }
  
  edge0 <- which(phy$edge[,2] == tip0)
  
  length0 <- phy$edge.length[edge0]
  
  ageofo0 <- mrcaage0 + length0
  
  # Create a new branch in which the branch length length1 is calculated so that fossil attaches in between its age and the time of origin of the parent lineage (controlled by edge.rel.length option).
  
  length1 <- (ageofo0-fossil.age)*edge.rel.length
  
  edge1 <- compute.brlen(stree(1, tip.label=fossil), length1)
  
  # Graft the branch into the tree:
  # The trick is to calculate the position so the fossil keeps its true age; which is the fossil age plus its branch length (minus mra age if the branch is not terminal)
  
  position <- fossil.age + length1 - mrcaage0
  
  new.tree <- bind.tree(phy, edge1, where=tip0, position=position)
  
  return(new.tree)
}

#function to return the terminal branch lengths of a tree
terminalLengths<-function(tree){
  tips<-tree$tip.label
  #http://blog.phytools.org/2016/02/extracting-terminal-edge-lengths-for.html
  ## first get the node numbers of the tips
  nodes<-sapply(tips,function(x,y) which(y==x),y=tree$tip.label)
  ## then get the edge lengths for those nodes
  edge.lengths<-setNames(tree$edge.length[sapply(nodes,
                                                 function(x,y) which(y==x),y=tree$edge[,2])],names(nodes))
  return(edge.lengths)
}


#cv function
cv <-  function(data) {
  
  return(sd(data) / mean(data) * 100)

  }
