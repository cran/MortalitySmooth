\name{plot.HMDdata}
\alias{plot.HMDdata}

\title{ Plot age and year specific demographic data }

\description{
  Plot one or two-dimensional mortality data from a \code{HMDdata}
  object created by \code{selectHMDdata}.  
}

\usage{
## S3 method for class 'HMDdata'
\method{plot}{HMDdata}(x, ...)
}

\arguments{
  \item{x}{ A \code{HMDdata} object. }
  \item{\dots}{ Other plotting parameters. Either passed to \code{plot}
  in a one-dimensional plot or passed to \code{levelplot} in a
  two-dimensional surface. } 
}

\details{
  The function needs \code{\link{HMDdata}} to be uploaded.
  
  A \code{HMDdata} object can be produced by \code{selectHMDdata}.
  
  The function recognizes the \code{HMDdata} object, its dimensions and
  its attributes (country, type of data and sex). Given the dimensions
  of the object the function plot either a simple unidimensional plot or
  a more complex shaded contour map over age and years.   
  
  Note that rates are automatically plotted in log-scale.
}

\value{
  None. Function produces a plot.
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{HMDdata}}, \code{\link{selectHMDdata}}. }

\examples{
## plot.HMDdata
data(HMDdata)
# 1-D plot
popDENf <- selectHMDdata("Denmark", "Population", "Females", 50:100, 2000) 
plot(popDENf, main="Danish female population, 2000", xlab="ages")
ratDENf <- selectHMDdata("Denmark", "Rates", "Females", 50:100, 2000)
plot(ratDENf, main="Danish female log death rates, 2000", xlab="ages")
# 2-D plot
years <- 1950:2000
ages <- 0:100
popJAPf <- selectHMDdata("Japan", "Population", "Females", ages, years)
plot(popJAPf, main="Japanese female population")
ratSWEf <- selectHMDdata("Sweden", "Rates", "Females")
plot(ratSWEf, main="Swedish female log death rates")
}

\keyword{hplot}
