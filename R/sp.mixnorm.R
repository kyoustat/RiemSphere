#' Mixture of Spherical Normal Distributions
#' 
#' @examples 
#' ## generate two-cluster data
#' mymu1 = c(0,0,0,1)  # center of class 1
#' mymu2 = c(-1,0,0,0) # center of class 2
#' mymu3 = c(0,1,0,0)  # center of class 3
#' 
#' x1 = rvmf(50, mymu1, kappa=10)
#' x2 = rvmf(50, mymu2, kappa=10)
#' x3 = rvmf(50, mymu3, kappa=10)
#' xx = rbind(x1,x2,x3)
#' 
#' ## apply clustering with different k values
#' sp.mixnorm(xx, k=2)
#' 
#' @export
sp.mixnorm <- function(x, k=2, n.start=5, maxiter=496, same.lambda=TRUE){
  ###################################################################
  ## Preprocessing 
  if (!check_datamat(x)){
    stop("* sp.mixnorm : input data should consist of at least two points on the unitsphere.")
  }

  myn = nrow(x)
  myp = ncol(x)
  myk = round(k)
  myiter = round(maxiter)
  
  if (myk >= (myn-1)){
    stop("* sp.mixnorm : the number of clusters is too big. Try with a smaller number.")
  }
  if (myk<2){
    stop("* sp.mixnorm : the number of clusters should be at least 2.")
  }
  
  ###################################################################
  ## Initialize all the parameters
  initlabel   <- (as.vector(stats::kmeans(x, k, nstart=n.start)$cluster)) # initialize the label
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
    loglkd.new = mixnorm.loglkd(x, par.mu, par.lambda, par.pi)
    
    
    # iteration over
    loglkd.old = loglkd.new
    print(paste("iteration ",it," with loglkd ",loglkd.old,sep=""))
  }
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
  
  output = 0
  for (k in 1:kk){
    output = output + sum(as.vector(dspnorm(x, as.vector(par.mu[k,]), lambda=par.lambda[k], log = TRUE)))
  }
  return(output)
}
