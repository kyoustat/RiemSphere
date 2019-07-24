#' Rayleigh Test of Uniformity
#' 
#' @export
testunif.Rayleigh <- function(x, method=c("Original","Modified")){
  ##############################################################
  # Preprocessing
  check_datamat(x)  # checking the datamatrix
  n = nrow(x)
  p = ncol(x)
  method = tolower(method)
  alldip = c("original","modified")
  method = match.arg(method, alldip)
  
  ##############################################################
  # Computation
  Rbar = aux_vmf_Rbar(x)
  if (aux_strcmp(method, "original")){
    hname   = "original"
    thestat = p*n*(Rbar^2)
    pvalue  = stats::pchisq(thestat, df=p, lower.tail=FALSE)
  } else if (aux_strcmp(method, "modified")){
    hname   = "modified"
    S       = p*n*(Rbar^2)
    thestat = (1-(1/(2*n)))*S + (1/(2*n*(p+2)))*(S^2)
    pvalue  = stats::pchisq(thestat, df=p, lower.tail=FALSE)
  } 
  
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  Ha      = "data is not uniformly distributed on sphere."
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}