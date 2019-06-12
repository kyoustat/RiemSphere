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
  # 3. x data
  check_datamat(x)
  
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
