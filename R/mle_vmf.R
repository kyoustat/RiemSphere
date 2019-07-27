#' MLE for von Mises-Fisher
#' 
#' 
#' @export
mle.vmf <- function(x, method=c("Banerjee","Song","Sra","Tanabe","Uniroot")){
  ########################################################################
  ## PREPROCESSING : check data
  check_datamat(x)
  if ((is.vector(x))||(nrow(x)==1)){
    stop("* mle.vmf : computing MLE requires more than one observations.")
  }
  method = tolower(method)
  alldip = c("banerjee","sra","tanabe","uniroot","song")
  method = match.arg(method, alldip)
  
  ########################################################################
  ## Compute 1 : extrinsic mean
  opt.mean = as.vector(colMeans(x))
  opt.mean = opt.mean/sqrt(sum(opt.mean^2))
  
  ########################################################################
  ## Compute 2 : kappa concentration
  opt.kappa <- switch (method,
    "banerjee" = vmf_2005banerjee(x),
    "sra"      = vmf_2012sra(x),
    "tanabe"   = vmf_2007tanabe(x),
    "uniroot"  = vmf_uniroot(x),
    "song"     = vmf_2012song(x)
  )
  
  return(opt.kappa)
}






######## SEVERAL KAPPA APPROXIMATION METHODS ##############################
# 1. 2005 Banerjee --------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_2005banerjee <- function(x){
  rbar <- aux_vmf_Rbar(x)
  p    <- ncol(x)
  
  term1 = rbar*(p-(rbar^2))
  term2 = 1-(rbar^2)
  return(term1/term2)
}


# 2. 2007 Tanabe ----------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_2007tanabe <- function(x){
  p       = ncol(x)
  rbar    = aux_vmf_Rbar(x)
  kap.old = vmf_2005banerjee(x)  # initialization
  
  maxiter = 123
  reltol  = 1e-7
  valtol  = 12345
  citer   = 1
  
  while (valtol > reltol){
    kap.new = rbar*kap.old/aux_vmf_Apk(p, kap.old)
    valtol  = abs(kap.new-kap.old)/kap.old
    citer = citer + 1
    
    kap.old = kap.new
    if (citer >= maxiter){
      break
    }
  }
  return(kap.old)
}


# 3. 2012 Sra -------------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_2012sra <- function(x){
  p    = ncol(x)
  rbar = aux_vmf_Rbar(x)
  kap0 = vmf_2005banerjee(x)     # initialization
  kap0Apk = aux_vmf_Apk(p, kap0) # A_p(kap0)
  
  # first iteration
  kap1 = kap0 - ((kap0Apk-rbar)/(1-(kap0Apk^2)-(((p-1)/kap0)*kap0Apk)))
  kap1Apk = aux_vmf_Apk(p, kap1)
  
  # second iteration
  kap2 = kap1 - ((kap1Apk-rbar)/(1-(kap1Apk^2)-(((p-1)/kap1)*kap1Apk)))
  return(kap2)
  
}


# 4. uniroot --------------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_uniroot <- function(x){
  rbar = aux_vmf_Rbar(x)
  p    = ncol(x)
  
  myfun <- function(kappa){
    return(aux_vmf_Apk(p, kappa)-rbar)
  }
  tgt = vmf_2005banerjee(x)
  return(stats::uniroot(myfun, lower=tgt/2, upper=tgt*2)$root)
}


# 5. 2012 Song ------------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_2012song <- function(x){
  p = ncol(x)
  kap0 = vmf_2005banerjee(x)
  rbar = aux_vmf_Rbar(x)
  
  # iteration 1.
  f0   = aux_vmf_Apk(p, kap0)-rbar
  df0  = aux_vmf_dApk(p, kap0)
  ddf0 = aux_vmf_d2Apk(p, kap0)
  kap1 = kap0 - (2*f0*df0)/(2*(df0^2) - (f0*ddf0))
  
  # iteration 2.
  f1   = aux_vmf_Apk(p, kap1)-rbar
  df1  = aux_vmf_dApk(p, kap1)
  ddf1 = aux_vmf_d2Apk(p, kap1)
  kap2 = kap1 - (2*f1*df1)/(2*(df1^2) - (f1*ddf1))
  return(kap2)
}



# ## simply test
# x = rvmf(50, c(1,rep(0,7)), 10)
# mle.vmf(x, method="banerjee")
# mle.vmf(x, method="song")
# mle.vmf(x, method="sra")
# mle.vmf(x, method="tanabe")
# mle.vmf(x, method="uniroot")
