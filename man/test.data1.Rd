\name{test.data1}
\alias{test.data1}
\docType{data}
\title{simulated two-component test data}

\description{simulated two-component test data used in function DeMixT}

\usage{test.data1, test.data1$y, test.data1$comp1, test.data1$truth}

\format{A list containing two matrices
  \describe{
    \item{\code{y}}{a matrix of expression data from mixed tumor samples}
    \item{\code{comp1}}{a matrix of expression data from reference component 1}
    \item{\code{truth}}{a matrix of true proportions, i.e., \eqn{\pi_1} and  \eqn{1 - \pi_1}}
  }
}

\examples{
data(test.data1)
test.data1
}
\keyword{datasets}