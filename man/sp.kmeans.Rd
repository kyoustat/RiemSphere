% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp.kmeans.R
\name{sp.kmeans}
\alias{sp.kmeans}
\title{k-Means Clustering for Data on Hypersphere}
\usage{
sp.kmeans(
  x,
  k = 2,
  init = c("kmeans", "random"),
  maxiter = 100,
  type = c("intrinsic", "extrinsic")
)
}
\arguments{
\item{x}{an \eqn{(n\times p)} row-stacked matrix for \eqn{\mathbb{S}^{p-1}}.}

\item{k}{the number of clusters to be found.}

\item{init}{method for setting initial configuration. \code{"random"} or \code{"kmeans"}.}

\item{maxiter}{maximum number of iterations to be run.}

\item{type}{type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.}
}
\value{
a named list containing
\describe{
\item{cluster}{length-\eqn{n} vector of class membership indices.}
}
}
\description{
k-Means Clustering for Data on Hypersphere
}
\examples{
## generate two-cluster data
mymu1 = c(0,0,0,1)  # center of class 1
mymu2 = c(-1,0,0,0) # center of class 2

x1 = rvmf(50, mymu1, kappa=10)
x2 = rvmf(50, mymu2, kappa=10)
xx = rbind(x1,x2)

## apply clustering with different k values
cl2 <- sp.kmeans(xx, k=2)
cl3 <- sp.kmeans(xx, k=3)
cl4 <- sp.kmeans(xx, k=4)

## compute 2-dimensional embedding for visualization
mds2d <- sp.mds(xx, ndim=2)
mdsx  <- mds2d$embed[,1]
mdsy  <- mds2d$embed[,2]

## compare via visualization
opar  <- par(mfrow=c(1,3), pty="s")
plot(mdsx, mdsy, col=cl2$cluster, main="k=2 means", pch=19)
plot(mdsx, mdsy, col=cl3$cluster, main="k=3 means", pch=19)
plot(mdsx, mdsy, col=cl4$cluster, main="k=4 means", pch=19)
par(opar)

}
