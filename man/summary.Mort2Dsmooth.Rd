\name{summary.Mort2Dsmooth}

\alias{summary.Mort2Dsmooth}
\alias{print.summary.Mort2Dsmooth}

\title{ Summary for Mort2Dsmooth objects }

\description{
  Summarizes the Poisson P-spline model fitted to a two-dimensional
  data. It returns various settings and measures. 
}
\usage{
\method{summary}{Mort2Dsmooth}(object, ...)
}

\arguments{
  \item{object}{ an object of class "Mort2Dsmooth", usually, a result of
  a call to \code{Mort2Dsmooth}. } 
  \item{\dots}{ further arguments passed to or from other methods. }
}

\details{
  \code{print.summary.Mort2Dsmooth} tries to be smart about formatting
  settings, outcomes, etc. After the matched call, the function presents
  several outcomes of the model, such as AIC, BIC, effective dimension,
  selected smoothing parameters, overdispersion parameter and a summary
  of the deviance residuals. The last lines show specifications and
  control parameters of the fitted model on both axes. 
}

\value{
  It produces an object of class \code{summary.Mort2Dsmooth} which
  contains exactly the same components of the associated
  \code{Mort2Dsmooth} object. 
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}}. }

\examples{
## selected data
ages <- 10:60
years <- 1950:2006
death <- selectHMDdata("Sweden", "Deaths", "Males",
                       ages = ages, years = years) 
exposure <- selectHMDdata("Sweden", "Exposures", "Males",
                          ages = ages, years = years)
## fit
fit <- Mort2Dsmooth(x=ages, y=years, Z=death,
                    offset=log(exposure),
                    method=3, lambdas=c(0.1, 1000))
## summary
summary(fit)
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
