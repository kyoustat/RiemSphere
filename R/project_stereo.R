#' Stereographic Projection
#' 
#' wikipedia : (1,0,0,0)
#' @export
project.stereo <- function(dat){
  ## PREPROCESSING : check data
  check_datamat(dat)
  
  ## APPLY THE FUNCTION AND RETURN
  return(apply(dat, 1, project.stereo.single))
}

#' Inverse of Stereographic Projection
#' 
#' @export
project.invstereo <- function(dat){
  return(apply(dat, 1, project.invstereo.single))
}


# single operations -------------------------------------------------------
#' @keywords internal
#' @noRd
project.stereo.single <- function(x){
  d = length(x)
  output = x[2:length(x)]/(1-x[1])
  return(output)
}
#' @keywords internal
#' @noRd
project.invstereo.single <- function(x){
  S2 = sum(x^2)
  return(c(((S2-1)/(S2+1)),((2*x)/(S2+1))))
}