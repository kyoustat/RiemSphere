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
  
  ## MAIN : use the ported version
  output = rvmf_port(n, mu, kappa)
  return(output)
}


#' @keywords internal
#' @noRd
rvmf_port <- function(n, mu, k){
  d <- length(mu)
  if (k > 0) {
    mu <- mu/sqrt(sum(mu^2))
    ini <- c(numeric(d - 1), 1)
    d1 <- d - 1
    v1 <- matrix(rnorm(n*d1), nrow=n)
    v <- v1 / sqrt( base::rowSums(v1^2) )
    b <- (-2 * k + sqrt(4 * k^2 + d1^2))/d1
    x0 <- (1 - b)/(1 + b)
    m <- 0.5 * d1
    ca <- k * x0 + (d - 1) * log(1 - x0^2)
    w <- rvmf_h(n,ca,d1,x0,m,k,b) # load directly from RcppArmadillo function
    S <- cbind(sqrt(1 - w^2) * v, w)
    if (sqrt(sum((ini-mu)^2)) < 100*.Machine$double.eps){
      A = diag(length(ini))
    } else {
      A <- aux_rotation(ini, mu)
    }
    x <- tcrossprod(S, A)
  }
  else {
    x = rvmf_uniform(n, mu)
  }
  return(x)
}


#' @keywords internal
#' @noRd
rvmf_uniform <- function(n, mu, k=0){
  d = length(mu)
  x1 = matrix(rnorm(n*d),nrow=n)
  x  = x1/sqrt(base::rowSums(x1^2))
  return(x)
}