#' Mixture of Spherical Normal Distributions
#' 
#' max loglkd; min ICs
#' 
#' @examples 
#' ## generate two-cluster data
#' mymu1 = c(0,0,0,1)        # center of class 1
#' mymu2 = c(-1,0,0,0)       # center of class 2
#' mymu3 = c(0,1/sqrt(2),-1/sqrt(2),0) # class 3
#' 
#' x1 = rvmf(50, mymu1, kappa=10)
#' x2 = rvmf(50, mymu2, kappa=15)
#' x3 = rvmf(50, mymu3, kappa=20)
#' xx = rbind(x1,x2,x3)
#' 
#' ## apply clustering with different k values
#' mix2 <- sp.mixnorm(xx, k=2)
#' mix3 <- sp.mixnorm(xx, k=3)
#' mix4 <- sp.mixnorm(xx, k=4)
#' 
#' ## compute 2-dimensional embedding for visualization
#' mds2d <- sp.mds(xx, ndim=2)
#' x  <- mds2d$embed[,1]
#' y  <- mds2d$embed[,2]
#' 
#' ## compare via visualization
#' opar <- par(mfrow=c(1,3),pty="s")
#' plot(x,y,col=mix2$cluster,main="k=2 mixture",pch=19)
#' plot(x,y,col=mix3$cluster,main="k=2 mixture",pch=19)
#' plot(x,y,col=mix4$cluster,main="k=2 mixture",pch=19)
#' par(opar)
#' 
#' ## extra visualization
#' \dontrun{
#' # run more models
#' mix5 <- sp.mixnorm(xx, k=5)
#' mix6 <- sp.mixnorm(xx, k=6)
#' mix7 <- sp.mixnorm(xx, k=7)
#' mix8 <- sp.mixnorm(xx, k=8)
#' 
#' # information criteria
#' xmat <- 2:8
#' ymat <- rbind(mix2$criteria, mix3$criteria, mix4$criteria, mix5$criteria, 
#'               mix6$criteria, mix7$criteria, mix8$criteria)
#' colnames(ymat) = colnames(mix2$criteria)
#' 
#' # plot with x11()
#' x11(width=12, height=6)
#' par(mfrow=c(2,4), pty="s")
#' plot(x, y, col=rainbow(8)[mix2$cluster], main="k=2 mixture", pch=19)
#' plot(x, y, col=rainbow(8)[mix3$cluster], main="k=3 mixture", pch=19)
#' plot(x, y, col=rainbow(8)[mix4$cluster], main="k=4 mixture", pch=19)
#' plot(x, y, col=rainbow(8)[mix5$cluster], main="k=5 mixture", pch=19)
#' plot(x, y, col=rainbow(8)[mix6$cluster], main="k=6 mixture", pch=19)
#' plot(x, y, col=rainbow(8)[mix7$cluster], main="k=7 mixture", pch=19)
#' plot(x, y, col=rainbow(8)[mix8$cluster], main="k=8 mixture", pch=19)
#' matplot(xmat, ymat, type="b", lwd=2, main="Info. Criteria",
#'         xlab="# clusters", col=1:4, lty=1:4, pch=1)
#' legend("topright", legend = colnames(ymat), col = 1:4, lty=1:4, pch=1)
#' } 
#' 
#' @export
sp.mixnorm <- function(x, k=2, init=c("kmeans","random"), maxiter=496, same.lambda=TRUE, 
                       version=c("soft","hard","stochastic")){
  ###################################################################
  ## Preprocessing 
  if (!check_datamat(x)){
    stop("* sp.mixnorm : input data should consist of at least two points on the unitsphere.")
  }

  myn = nrow(x)
  myp = ncol(x)
  myk = round(k)
  myiter = round(maxiter)
  myinit = match.arg(init)
  myvers = match.arg(version)
  
  if (myk >= (myn-1)){
    stop("* sp.mixnorm : the number of clusters is too big. Try with a smaller number.")
  }
  if (myk<2){
    stop("* sp.mixnorm : the number of clusters should be at least 2.")
  }
  
  ###################################################################
  ## Initialize all the parameters
  if (all(myinit=="random")){
    initlabel = sample(c(1:myk, sample(1:myk, myn-myk, replace = TRUE)))
  } else {
    initlabel = stats::kmeans(x, myk, nstart=round(5))$cluster # label  
  }
  par.eta <- array(0,c(myn,myk))
  for (i in 1:myn){
    par.eta[i,initlabel[i]] = 1
  }
  par.mu = mixnorm.frechet(x, par.eta)
  d2mat  = (sp.pdist2.internal(x, par.mu, type="intrinsic")^2)
  if (same.lambda){
    par.lambda = mixnorm.lambda.homo(d2mat, par.eta, myp)
  } else {
    par.lambda = mixnorm.lambda.hetero(d2mat, par.eta, myp)
  }
  par.pi = base::colSums(par.eta)/myn
  loglkd.old = mixnorm.loglkd(x, par.mu, par.lambda, par.pi)

  ###################################################################
  ## Let's Iterate
  for (it in 1:myiter){
    # E-Step
    par.eta = mixnorm.eta(x, par.mu, par.lambda, par.pi)
    
    # H/S-Step by Option
    if (aux_strcmp(myvers, "hard")){
      par.eta = mixnorm.hard(par.eta)
    } else if (aux_strcmp(myvers, "stochastic")){
      par.eta = mixnorm.stochastic(par.eta)
    }
    
    
    # M-Step Parameter Update
    #   Preliminary. d2mat
    d2mat = (sp.pdist2.internal(x, par.mu, type="intrinsic")^2)
    #   M1. lambda / concentration
    if (same.lambda){
      par.lambda = mixnorm.lambda.homo(d2mat, par.eta, myp)
    } else {
      par.lambda = mixnorm.lambda.hetero(d2mat, par.eta, myp)
    }
    #   M2. proportion
    par.pi = base::colSums(par.eta)/myn
    #   M3. mu / centers
    par.mu = mixnorm.frechet(x, par.eta)
    #   Mextra for computing log-likelihood
    loglkd.new  = mixnorm.loglkd(x, par.mu, par.lambda, par.pi)
    loglkd.diff = loglkd.new - loglkd.old
    # iteration over
    loglkd.old = loglkd.new
    if ((it > 5)&&(loglkd.diff <= 0)){
      break
    }
  }
  
  ###################################################################
  # Information Criterion
  if (!same.lambda){
    par.k = ((myp-1)*myk) + myk + (myk-1)  
  } else {
    par.k = ((myp-1)*myk) + 1   + (myk-1)  
  }
  AIC = -2*loglkd.old + 2*par.k
  BIC = -2*loglkd.old + par.k*log(myn)
  HQIC = -2*loglkd.old + 2*par.k*log(log(myn))
  AICc = AIC + (2*(par.k^2) + 2*par.k)/(myn-par.k-1)
  
  infov = matrix(c(AIC, AICc, BIC, HQIC), nrow=1)
  colnames(infov) = c("AIC","AICc","BIC","HQIC")
  
  
  ###################################################################
  # Return
  output = list()
  output$cluster  = base::apply(par.eta, 1, aux_whichmax)
  output$loglkd   = loglkd.old  # max loglkd
  output$criteria = infov       # min AIC/AICc/BIC/HQIC
  output$parameters = list(pro=par.pi, centers=par.mu, concentration=par.lambda)
  output$membership = par.eta
  return(output)
}


#   -----------------------------------------------------------------------
#' 1. compute weighted frechet mean
#' @keywords internal
#' @noRd
mixnorm.frechet <- function(dat, eta){
  n = nrow(dat)
  p = ncol(dat)
  k = ncol(eta)
  output = array(0,c(k,p))

    for (kk in 1:k){
    if (sum(eta[,kk]==1)==1){
      output[kk,] = dat[which(eta[,kk]==1),]
    } else {
      output[kk,] = as.vector(aux_wfrechet(dat, weight=as.vector(eta[,kk])))  
    }
  }
  return(output)
}
#' 2. compute lambda for homogeneous/heterogeneous
#' @keywords internal
#' @noRd
mixnorm.lambda.homo <- function(d2mat, eta, myp){
  k = ncol(eta)
  term1 = base::sum(d2mat*eta)
  term2 = base::sum(eta)
  
  myfun <- function(lambda){ # function to be minimized
    return((lambda*term1/2) + (term2*(log(dspnorm.constant(lambda,myp)))))
  }
  
  output = stats::optimize(myfun, interval=c(10*.Machine$double.eps,12345), maximum = FALSE)$minimum
  return(rep(output, k))
}
#' @keywords internal
#' @noRd
mixnorm.lambda.hetero <- function(d2mat, eta, myp){
  kk = ncol(eta)
  outvec = rep(0,kk)
  d2meta = d2mat*eta
  
  for (k in 1:kk){
    term1 = as.vector(sum(d2meta[,k]))
    term2 = as.vector(sum(eta[,k]))
    
    myfun <- function(lambda){
      return(((lambda*term1/2) + (term2*(log(dspnorm.constant(lambda,myp))))))
    }
    
    outvec[k] = stats::optimize(myfun, interval=c(10*.Machine$double.eps,12345), maximum = FALSE)$minimum
  }
  return(outvec)
}
#' 3. compute eta parameter
#' @keywords internal
#' @noRd
mixnorm.eta <- function(x, par.mu, par.lambda, par.pi){
  kk = length(par.pi)
  nn = nrow(x)
  
  par.eta = array(0,c(nn,kk))
  for (k in 1:kk){
    par.eta[,k] = as.vector(dspnorm(x, as.vector(par.mu[k,]), lambda=par.lambda[k], log = FALSE))
  }
  output = par.eta/base::rowSums(par.eta)
  return(output)
}
#' 4. compute log-likelihood
#' @keywords internal
#' @noRd
mixnorm.loglkd <- function(x, par.mu, par.lambda, par.pi){
  kk = length(par.pi)
  nn = nrow(x)
  
  evals = array(0,c(nn,kk))
  for (k in 1:kk){
    evals[,k] = as.vector(dspnorm(x, as.vector(par.mu[k,]), par.lambda[k], log = FALSE))*par.pi[k]
  }
  return(base::sum(base::log(base::rowSums(evals)))) # following my current note, works fine !
}
#' 5. version 'hard' 
#' @keywords internal
#' @noRd
mixnorm.hard <- function(eta){
  n = nrow(eta)
  p = ncol(eta)
  idmax  = base::apply(eta, 1, aux_whichmax)
  output = array(0,c(n,p))
  for (i in 1:n){
    output[i,idmax[i]] = 1
  }
  return(output)
}
#' 6. version 'stochastic' 
#' @keywords internal
#' @noRd
mixnorm.stochastic <- function(eta){
  n = nrow(eta)
  k = ncol(eta)
  vec1k = (1:k)
  output = array(0,c(n,k))
  for (i in 1:n){
    output[i,base::sample(vec1k, 1, prob=as.vector(eta[i,]))] = 1
  }
  return(output)
}



# rm(list=ls())
# library(RiemSphere)
# data("spdat.orbital")
# 
# dmat = sp.pdist(spdat.orbital)
# mds2 = DAS::cmds(dmat)$embed
# x = as.vector(mds2[,1])
# y = as.vector(mds2[,2])
# 
# c2 = sp.kmeans(spdat.orbital, k=2)
# c3 = sp.kmeans(spdat.orbital, k=3)
# c4 = sp.kmeans(spdat.orbital, k=4)
# c5 = sp.kmeans(spdat.orbital, k=5)
# 
# graphics.off()
# x11()
# par(mfrow=c(2,2),pty="s")
# plot(x,y,col=c2$cluster,pch=19,main="k=2")
# plot(x,y,col=c3$cluster,pch=19,main="k=3")
# plot(x,y,col=c4$cluster,pch=19,main="k=4")
# plot(x,y,col=c5$cluster,pch=19,main="k=5")