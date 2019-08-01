#' One-sample Location Test with Log-likelihood Ratio Test
#' 
#' @export
spmean1.lrt <- function(x, mu0, kappa0=NULL){
  ##############################################################
  # Preprocessing
  check_datamat(x)  # checking the datamatrix
  check_datamat(matrix(mu0, nrow=1))
  n = nrow(x)
  p = ncol(x)
  if (length(mu0)!=p){
    stop("* spmean1.lrt : 'mu' should have length of 'ncol(x)'.")
  }
  
  ##############################################################
  # BRANCHING
  if ((length(kappa0)==1)&&(!is.null(kappa))){      # Case1 : Known
    hname   = "Known"
    Xbar = as.vector(colMeans(x))
    Rbar = sqrt(sum(Xbar^2))
    Cbar = sum(Xbar*mu0)
    thestat =  2*n*kappa0*(Rbar - Cbar)
    pvalue  = stats::pchisq(thestat, df=p, lower.tail=FALSE)
  } else {                                          # Case2 : Unknown
    hname   = "Unknown"
    Xbar = as.vector(colMeans(x))
    Rbar = sqrt(sum(Xbar^2))
    ktilde = aux_vmf_ReML_kappa(x, mu0)
    khat   = mle.vmf(x, "nmle1")$kappa
    thestat = n*(khat*Rbar - ktilde*sum(mu0*Xbar) - aux_vmf_apk(p,khat) + aux_vmf_apk(p,ktilde))
    pvalue  = stats::pchisq(thestat, df=p, lower.tail=FALSE)
  }
  
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  Ha      = "true mean is different from mu0."
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}

# # example
# mymu = rnorm(4)
# mymu = mymu/sqrt(sum(mymu^2))
# myx  = RiemSphere::rvmf(100, mymu, kappa=1)
# spmean1.lrt(myx, mymu, kappa0=NULL)
# spmean1.lrt(myx, mymu, kappa0=1)

