#' The von Mises-Fisher Distribution
#' 
#' Group of functions Description section
#' 
#' Group of functions Details paragraph.
#'
#' @section After Arguments and Value sections:
#' Despite its location, this actually comes after the Arguments and Value sections.
#' Also, don't need to use null, could annotate first function, and then
#' using function name as the groupBy name is more intuitive.
#' 
#' @param x a param for toBar and notToBar
#' @param y a param just for notToBar
#' @return Hard to have one return section for all functions,
#' might want to have a manual list here.
#' @name vmf
NULL

#' @rdname vmf
#' @export
dvmf <- function(x, mu, kappa=1, log = FALSE){
  ## PREPROCESSING
  # 1. mu : unit vector
  check_unitvec(mu)
  p = length(mu)
  # 2. lambda 
  check_num_nonneg(kappa)
  
  ## EVALUATION
  #   1. normalizing term
  nterm = (kappa^((p/2)-1))/(((2*pi)^(p/2))*besselI(kappa,((p/2)-1)))
  #   2. case branching
  if (is.vector(x)){
    output = exp(kappa*sum(mu*x))
  } else {
    output = exp(kappa*as.vector(x%*%mu))
  }
  #   3. scale by normalizing constant and RETURN
  if (log){
    return(log(output)+log(nterm))
  } else {
    return(output*nterm)
  }
}
# 
# # compare with vmf.density
# niter = 12345
# val1 = rep(0,niter)
# val2 = rep(0,niter)
# mydat = array(0,c(niter,3))
# mymu = rnorm(3)
# mymu = mymu/sqrt(sum(mymu^2))
# mykappa = abs(rnorm(1))
# for (i in 1:niter){
#   y = rnorm(3)
#   y = y/sqrt(sum(y^2))
#   mydat[i,] = y
#   val1[i] = vmf.density(y, mykappa, mymu)
#   val2[i] = dvmf(y, mymu, mykappa)
# }
# val3 = dvmf(mydat, mymu, mykappa)


#   -----------------------------------------------------------------------
#' @rdname vmf
#' @param n the number
#'
#' Could put here, but easiest if all params are
#' described in the same place, with page documentation, so none get repeated.
#' @export
rvmf <- function(n, mu, kappa=1){
  ## PREPROCESSING
  # 1. mu : unit vector
  check_unitvec(mu)
  # 2. lambda 
  check_num_nonneg(kappa)
  
  ## MAIN : use Rfast package
  output = Rfast::rvmf(n, mu, kappa)
  return(output)
}
