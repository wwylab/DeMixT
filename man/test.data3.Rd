\name{test.data3}
\alias{test.data3}
\docType{data}
\title{
three-component mixed cell line test data
}

\description{
three-component mixed cell line test data used in function DeMixT
}

\usage{test.data3, test.data3$y, test.data3$comp1, test.data3$comp2}
\format{
  A list containing three matrices
  \describe{
    \item{\code{y}}{a matrix of expression data from mixed tumor samples}
    \item{\code{comp1}}{a matrix of expression data from reference component 1}
    \item{\code{comp2}}{a matrix of expression data from reference component 2}
  }
}
\examples{
data(test.data3)
test.data3
}
\keyword{datasets}