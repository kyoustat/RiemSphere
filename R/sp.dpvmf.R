#' Dirichet Process Mixture of von Mises-Fisher Distributions
#' 
#' @param data an \eqn{(n\times d)} data matrix where \eqn{x_i \in \mathcal{S}^{d-1}}.
#' @param kappa0 d
#' @param mu0 d
#' @param a d
#' @param b d
#' @param w0 d
#' @param mh.kappa stepsize for MH algorithm.
#' @param maxiter d
#' @param print.progress d
#' 
#' @return a length-\code{maxiter} list whose elements are lists containing \describe{
#' \item{a}{b}
#' }
#' 
#' @export
sp.dpvmf <- function(data, kappa0){
  ############################################################
  # Preprocessing
  if (!check_datamat(data)){
    stop("* sp.dpvmf : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
}