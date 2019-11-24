#' Normal Vectors to the Orbital Planes of the 9 Planets
#' 
#' \emph{Lena} is probably one of the most well-known example in image processing and computer vision.
#' Well, here is a brief introduction on \href{http://www.lenna.org/}{the story of Lena}.
#'
#'
#' @docType data
#' @usage data(data.orbital)
#' @format matrix of size \eqn{(9\times 3)}
#' @keywords datasets
#' @references Gonzalez, Rafael C. and Woods, Richard E. (2017) \emph{Digital Image Processing} (4th ed.). ISBN 0133356728.
#'
#' @source \href{http://sipi.usc.edu/database/?volume=misc}{USC SIPI Image Database}
#'
#' @examples
#' ## load the data
#' data(spdat.orbital)
#' 
#' ## visualize using multidimensional scaling
#' dmat = sp.pdist(spdat.orbital)
#' mds2 = DAS::cmds(dmat, ndim=2)
#' 
#' x = mds2$embed[,1]
#' y = mds2$embed[,2]
#' 
#' plot(x, y, pch=19, main="orbital normal vectors",
#'      xlim=c(min(x)-0.05,max(x)+0.05), ylim=c(min(y)-0.05,max(y)+0.05))
#' text(x+0.01, y+0.01, 1:9)
"spdat.orbital"