\name{summary.Mort1Dsmooth}

\alias{summary.Mort1Dsmooth}
\alias{print.summary.Mort1Dsmooth}

\title{ Summary for Mort1Dsmooth objects }

\description{
  Summarizes the Poisson P-spline model fitted to a unidimensional
  data. It returns various settings and measures. 
}
\usage{
\method{summary}{Mort1Dsmooth}(object, ...)
}

\arguments{
  \item{object}{ an object of class "Mort1Dsmooth", usually, a result of
  a call to \code{Mort1Dsmooth}. } 
  \item{\dots}{ further arguments passed to or from other methods. }
}

\details{
  \code{print.summary.Mort1Dsmooth} tries to be smart about formatting
  settings, outcomes, etc. After the matched call, the function presents
  several outcomes of the model, such as AIC, BIC, effective dimension,
  selected smoothing parameter, overdispersion parameter and a summary
  of the deviance residuals. The last lines show specifications and
  control parameters of the fitted model. 
}

\value{
  It produces an object of class \code{summary.Mort1Dsmooth} which
  contains exactly the same components of the associated
  \code{Mort1Dsmooth} object. 
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort1Dsmooth}}. } 

\examples{
## selected data
years <- 1970:2006
death <- selectHMDdata("Sweden", "Deaths", "Females",
                       ages = 0, years = years)
exposure <- selectHMDdata("Sweden", "Exposures", "Females",
                          ages = 0, years = years)
## fit
fit <- Mort1Dsmooth(x=years, y=death, offset=log(exposure),
                    method=3, lambda=30)
## summary
summary(fit)
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
