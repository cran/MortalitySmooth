\name{plot.Mort2Dsmooth}
\alias{plot.Mort2Dsmooth}

\title{ Plot of the Mort2Dsmooth outcomes }

\description{
  It produces two level plots, side-by-side, where the x- and y-axis are
  the predictors (\code{x} and \code{y}) of the given
  \code{Mort2Dsmooth} object and the numeric response is either the
  matrices of actual and fitted log-rates, or the matrices of actual and
  fitted counts. 
}

\usage{
## S3 method for class 'Mort2Dsmooth'
\method{plot}{Mort2Dsmooth}(x, type = c("logrates", "deaths"), ...)
}

\arguments{
  \item{x}{ An object of class "Mort2Dsmooth", usually, a result of
  a call to \code{Mort2Dsmooth}. } 
  \item{type}{ The type of plot which should be returned. The
    alternatives are: \code{logrates} (default) and \code{death}. }
  \item{\dots}{ Other plotting parameters passed to \code{levelplot}. }
}

\details{
  The function provides a simple tool for checking the outcomes of a
  \code{Mort2Dsmooth} object. In a Poisson setting, one would commonly
  look at rates in log-scale, but actual and fitted counts can be
  plotted too. 
  
  The function uses \code{\link{levelplot}} from package \code{lattice} to construct the
  two level plots side-by-side. Plotting parameters can not be
  changed. A rainbow palette with 9 colors is chosen so that lower
  values are red, followed by orange, yellow, green, blue and indigo
  with the highest values plotted in violet. The breaks for the level
  plot are chosen according to the deciles of distributions of actual
  and fitted log-rates (counts). 
}

\value{
  None. Function produces a plot.
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}} for computing
  \code{Mort2Dsmooth.object}. } 

\examples{
# selected data
ages <- 50:100
years <- 1970:2006
death <- selectHMDdata("Denmark", "Deaths", "Females",
                       ages = ages, years = years) 
exposure <- selectHMDdata("Denmark", "Exposures", "Females",
                          ages = ages, years = years)
# fit
fit <- Mort2Dsmooth(x=ages, y=years, Z=death, offset=log(exposure),
                    method=3, lambdas=c(100,500)) 
# plotting log-death rates (default)
plot(fit)
# plotting death counts
plot(fit, type="death")
}

\keyword{hplot}
\keyword{regression}
\keyword{smooth}
