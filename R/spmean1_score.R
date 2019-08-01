#' One-sample Location Test with Score Test
#' 
#' @export
spmean1.score <- function(x, mu0, kappa0=NULL){
  ##############################################################
  # Preprocessing
  check_datamat(x)  # checking the datamatrix
  check_datamat(matrix(mu0, nrow=1))
  n = nrow(x)
  p = ncol(x)
  if (length(mu0)!=p){
    stop("* spmean1.score : 'mu' should have length of 'ncol(x)'.")
  }
  
  ##############################################################
  # BRANCHING
  if ((length(kappa0)==1)&&(!is.null(kappa))){      # Case1 : Known
    hname   = "Known Score Test"
    Xbar    = as.vector(colMeans(x))
    thestat = n*kappa0*(sum(as.vector((diag(p)-outer(mu0,mu0))%*%Xbar)^2))/aux_vmf_Apk(p,kappa0)
    pvalue  = stats::pchisq(thestat, df=p-1, lower.tail=FALSE)
  } else {                                          # Case2 : Unknown
    hname   = "Unknown Score Test"
    kappa0  = aux_vmf_ReML_kappa(x, mu0)
    Xbar    = as.vector(colMeans(x))
    thestat = n*kappa0*(sum(as.vector((diag(p)-outer(mu0,mu0))%*%Xbar)^2))/aux_vmf_Apk(p,kappa0)
    pvalue  = stats::pchisq(thestat, df=p-1, lower.tail=FALSE)
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
# spmean1.score(myx, mymu, kappa0=NULL)
# spmean1.score(myx, mymu, kappa0=1)
