\name{tpe1d}
\alias{tpe1d}
\title{One-dimensional tree preserving embedding}
\description{Finds a one-dimensional tree preserving embedding of a dissimilarity matrix.}
\usage{tpe1d(d, verbose=TRUE)}
\arguments{
\item{d}{A \code{\link[stats]{dist}} object with the dissimilarity matrix.}
\item{verbose}{A logical value indicating whether to print progress.}
}
\details{Implements the greedy approximation for one-dimensional tree preserving embedding.}
\value{A matrix with the embedding coordinates.}
\references{Shieh AD, Hashimoto TB, Airoldi EM (2011) Tree preserving embedding. \emph{Proc Natl Acad Sci USA} 108:16916--16921.}
