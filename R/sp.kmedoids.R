#' k-Medoids Clustering for Data on Hypersphere
#' 
#' 
#' @export
sp.kmedoids <- function(x, k=2, type=c("intrinsic","extrinsic")){
  ############################################################
  # Preprocessing
  if (!check_datamat(x)){
    stop("* sp.kmedoids : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
  myk    = round(k)
  mytype = match.arg(type)
  
  ############################################################
  # Compute Pairwise Distances and do PAM
  dmat   = sp.pdist.internal(x, type=mytype, as.dist=TRUE)
  output = RiemBaseExt::rclust.kmedoids(dmat, k=myk)
  return(output)
}