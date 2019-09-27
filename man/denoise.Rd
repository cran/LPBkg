\name{denoise}
\alias{denoise}
\title{Coefficients of the denoised comparison density estimator}
\description{
Selects the largest coefficients according to the AIC or BIC criterion.
}

\usage{
denoise(LP, n, method)
}
\arguments{
\item{LP}{Original vector of coefficients estimates. See details.}

\item{n}{The dimension of the sample on which the estimates in \code{LP} have been obtained.}

\item{method}{Either ``AIC'' or ``BIC''. See details.}
}
\details{Give a vector of \code{M} coefficient estimates, the largest is selected according to the AIC or BIC criterion as described in Algeri, 2019 and Mukhopadhyay, 2017. }

\value{
Selected coefficient estimates.
}

\references{
S. Algeri, 2019. Detecting new signals under background mismodelling.	<arXiv:1906.06615>.

S. Mukhopadhyay, 2017. Large-scale mode identification and data-driven sciences. Electronic Journal of Statistics 11 (2017), no. 1, 215--240.
}
\author{Sara Algeri}


\seealso{\code{\link{Legj}}.}

\examples{
#generating data
x<-rnorm(1000,10,7)
xx<-x[x>=10 & x<=20]

#create suitable postulated quantile function
G<-pnorm(20,5,15)-pnorm(10,5,15)
g<-function(x){dnorm(x,5,15)/G}

#Vectorize quantile function
g<-Vectorize(g)
u<-g(xx)

Mmax=20
S<- as.matrix(Legj(u=u,m=Mmax))
n<-length(u)

LP <- apply(S,FUN="mean",2)

denoise(LP,n=n,method="AIC")
}

\keyword{Denoised estimator.}
