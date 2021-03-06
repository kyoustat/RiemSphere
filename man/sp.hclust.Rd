% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp.hclust.R
\name{sp.hclust}
\alias{sp.hclust}
\title{Hierarchical Clustering for Data on Hypersphere}
\usage{
sp.hclust(
  x,
  type = c("intrinsic", "extrinsic"),
  method = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2",
    "centroid", "median"),
  members = NULL
)
}
\arguments{
\item{x}{an \eqn{(n\times p)} row-stacked matrix for \eqn{\mathbb{S}^{p-1}}.}

\item{type}{type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.}

\item{method}{the agglomeration method to be used. This must be (an unambiguous abbreviation of) one of \code{"single"},
\code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.}

\item{members}{\code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.}
}
\value{
an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details.
}
\description{
Hierarchical Clustering for Data on Hypersphere
}
\examples{
## generate two-cluster data
mymu1 = c(0,0,0,1)  # center of class 1
mymu2 = c(-1,0,0,0) # center of class 2

x1 = rvmf(20, mymu1, kappa=5)
x2 = rvmf(20, mymu2, kappa=5)
xx = rbind(x1,x2)

## apply hierarchical clustering with different methods
hc1 <- sp.hclust(xx, method="single")
hc2 <- sp.hclust(xx, method="complete")
hc3 <- sp.hclust(xx, method="average")

## visualize
\dontrun{
opar <- par(mfrow=c(1,3), pty="s")
plot(hc1, main="'single'")
plot(hc2, main="'complete'")
plot(hc3, main="'average'")
par(opar)
}

}
