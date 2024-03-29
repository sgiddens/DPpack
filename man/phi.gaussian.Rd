% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatModels.R
\name{phi.gaussian}
\alias{phi.gaussian}
\title{Transform Function for Gaussian Kernel Approximation}
\usage{
phi.gaussian(x, theta)
}
\arguments{
\item{x}{Vector or matrix corresponding to one row of the dataset X.}

\item{theta}{Randomly sampled prefilter vector of length n+1, where n is the
length of x.}
}
\value{
Mapped value corresponding to one element of the transformed space.
}
\description{
This function maps an input data row x with a given prefilter to an output
value in such a way as to approximate the Gaussian kernel
\insertCite{chaudhuri2011}{DPpack}.
}
\examples{
  x <- c(1,2,3)
  theta <- c(0.1, 1.1, -0.8, 3)
  phi.gaussian(x, theta)

}
\references{
\insertRef{chaudhuri2011}{DPpack}
}
\keyword{internal}
