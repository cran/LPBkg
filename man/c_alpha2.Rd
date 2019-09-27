\name{c_alpha2}
\alias{c_alpha2}
\title{Approximated quantiles}
\description{
Approximates the quantiles of the supremum of the comparison density estimator using tube formulae and assuming that $H_0$ is true.
}
\usage{c_alpha2(M, IDs, alpha = 0.05, c_interval = c(1, 10))}

\arguments{
\item{M}{The size of the polynomial basis used to estimate the comparison density.}

\item{IDs}{The IDs of the polynomial terms to be used out of the \code{M} considered.}

\item{alpha}{Desired significance level.}

\item{c_interval}{Lower and  upper bounds for the quantile being computed.}
}
\value{
Approximated quantile of order \code{1-alpha} of the supremum of the comparison density estimator.
}
\references{S. Algeri, 2019. Detecting new signals under background mismodelling	<arXiv:1906.06615>.

L.A. Wasserman, 2005. All of Nonparametric Statistics. Springer Texts in Statistics.}

\author{Sara Algeri}

\seealso{\code{\link{dhatL2}}.}
\examples{
c_alpha2(5, c(2,4), alpha = 0.05, c_interval = c(1, 10))
}



\keyword{ Tube formulae }

