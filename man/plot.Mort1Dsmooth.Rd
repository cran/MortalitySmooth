\name{plot.Mort1Dsmooth}
\alias{plot.Mort1Dsmooth}

\title{ Plot of the Mort1Dsmooth outcomes }

\description{
  It produces a plot where the x-axis is the predictor of the given
  \code{Mort1Dsmooth} object and the y-axis is both actual and
  fitted log-rates, or actual and fitted counts. 
}

\usage{
\method{plot}{Mort1Dsmooth}(x, Type = c("logrates", "deaths"), ...)
}

\arguments{
  \item{x}{ an object of class "Mort1Dsmooth", usually, a result of
  a call to \code{Mort1Dsmooth}. } 
  \item{Type}{ the type of plot which should be returned. The
    alternatives are: \code{logrates} (default) and \code{death}. }
  \item{\dots}{ other plotting parameters passed to \code{plot}. Not in used. }
}

\details{
  The function provides a simple tool for checking the outcomes of a
  \code{Mort1Dsmooth} object. In a Poisson setting, one would commonly
  look at rates in log-scale, but actual and fitted counts can be
  plotted too. 

  Plotting parameters are given and cannot be changed. 
}

\value{
  None. Function produces a plot
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort1Dsmooth}} for computing
  \code{Mort1Dsmooth.object}. }  

\examples{
## selected data
years <- 1930:2006
death <- selectHMDdata("Denmark", "Deaths", "Females",
                       ages = 60, years = years)
exposure <- selectHMDdata("Denmark", "Exposures", "Females",
                          ages = 60, years = years)
## fit
fit <- Mort1Dsmooth(x=years, y=death, offset=log(exposure),
                    method=3, lambda=100)
## plotting actual and fitted data
par(mfrow=c(1,2))
plot(fit)
plot(fit, "deaths")
par(mfrow=c(1,1))
}

\keyword{hplot}
\keyword{regression}
\keyword{smooth}
