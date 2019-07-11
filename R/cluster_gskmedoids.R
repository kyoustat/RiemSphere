#' Clustering with Geodesic Spherical k-Medoids Algorithm
#' 
#' 
#' @export
gskmedoids <- function(x, k=2, maxiter=100){
  ## Preprocessing
  check_datamat(x)  # checking the datamatrix
  if (is.vector(x)||(nrow(x)==1)){
    stop("* gskmedoids : input data should consist of at least two points on the unitsphere.")
  }
  n = nrow(x)
  if (k >= (n-1)){
    stop("* gskmedoids : the number of clusters is too big. Try with a smaller number.")
  }
  if (k<2){
    stop("* gskmeans : the number of clusters should be at least 2.")
  }
  
  ## Compute Pairwise Distance
  pdmat <- stats::as.dist(aux_dist_MtoM(x))
  outee <- cluster::pam(pdmat, k)
  
  ## Report the Results
  output = list()
  output$cluster = outee$clustering
  output$centers = x[outee$id.med,]
  return(output)
}