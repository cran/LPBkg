\name{BestM}
\alias{BestM}
\title{Chooses the size of the polynomial basis}
\description{
Computes the deviance p-values considering different sizes of the polynomial basis and selects the one for which the deviance p-value is the smallest.
}
\usage{BestM(data, g, Mmax = 20, range = c(min(data),max(data)))}

\arguments{
\item{data}{A vector of data. See details. }

\item{g}{The postulated model from which we want to assess if deviations occur.}

\item{Mmax}{The maximum size of the polynomial basis from which a suitable value \code{M} is selected (the default is 20). See details.}

\item{range}{Range of the data/ search region considered. }
}
\details{
The argument \code{data} collects the data for which we want to test if deviations occur from the postulated model specified in the argument \code{g}. As in Algeri, 2019, the sample specified under \code{data} corresponds to the source-free sample in the background calibration phase and to the physics sample in the signal search phase.
The value \code{M} selected determines the smoothness of the estimated comparison density, with smaller values of M leading to smoother estimates. The deviance test is used to select the value \code{M} which leads to the most significant deviation from the postulated model. The default value for \code{Mmax} is set to 20. Notice that numerical issues may
arise for larger values of \code{Mmax}.}
\value{
\item{pvals}{The deviance test p-value obtained for each values of \code{M} (from 1 to \code{Mmax}) considered.}
\item{minp}{The minimum value of the deviance p-values observed.}
\item{Msel}{The value of \code{M} at which the minimum deviance p-values is achieved.}
}
\references{S. Algeri, 2019. Detecting new signals under background mismodelling	<arXiv:1906.06615>. }
\author{Sara Algeri}
\seealso{
\code{\link{denoise}}, \code{\link{dhatL2}}.
}

\examples{
#Generating data
x<-rnorm(1000,10,7)
data<-x[x>=10 & x<=20]

#Create suitable postulated quantile function of data
G<-pnorm(20,5,15)-pnorm(10,5,15)
g<-function(x){dnorm(x,5,15)/G}

Mmax=10
range=c(10,20)

BestM(data,g,Mmax,range)

}

\keyword{ Truncate legendre series }
\keyword{ Model selection }
