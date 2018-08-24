\name{test.data2}
\alias{test.data2}
\docType{data}
\title{simulated three-component test data}

\description{simulated three-component test data used in function DeMixT}

\usage{test.data2, test.data2$y, test.data2$comp1, test.data2$comp2, test.data2$truth}

\format{A list containing three matrices
  \describe{
    \item{\code{y}}{a matrix of expression data from mixed tumor samples}
    \item{\code{comp1}}{a matrix of expression data from reference component 1}
    \item{\code{comp2}}{a matrix of expression data from reference component 2}
    \item{\code{truth}}{a matrix of true proportions, i.e., \eqn{\pi_1}, \eqn{\pi_2}, and  \eqn{1 - \pi_1 - \pi_2}}
  }
}

\examples{
data(test.data2)
test.data2
}

\keyword{datasets}