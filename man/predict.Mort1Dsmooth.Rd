\name{predict.Mort1Dsmooth}
\alias{predict.Mort1Dsmooth}

\title{ Predict Method for 1D P-splines Fits }

\description{
  Obtains predictions and optionally estimates standard errors of those
  predictions from a fitted \code{Mort1Dsmooth} object. 
}

\usage{
## S3 method for class 'Mort1Dsmooth'
\method{predict}{Mort1Dsmooth}(object, newdata = NULL, 
        type = c("link", "response"), 
        se.fit = FALSE, ...)
}

\arguments{
  \item{object}{ An object of class "Mort1Dsmooth", usually, a result of
  a call to \code{Mort1Dsmooth}. } 
  \item{newdata}{ Optionally, a vector in which to look for x with which
  to predict. If omitted, the fitted linear predictors are used. } 
  \item{type}{ The type of prediction required. The default ("link") is
  on the scale of the linear predictors; the alternative "response" is
  on the scale of the response variable. } 
  \item{se.fit}{ Logical switch indicating if standard errors are
    required. Default: \code{FALSE}. }
  \item{\dots}{ Other predict parameters to passed to \code{predict}.} 
}

\details{
  If \code{newdata} is omitted the predictions are based on the data
  used for the fit. Note that, in common with other prediction
  functions, any offset supplied as an argument is always ignored when
  predicting, unlike offsets specified in modelling. 
  
  Forecast is not possible, therefore \code{newdata} has to be include
  within the range of the original \code{x} and \code{y}. 
}

\value{
  If \code{se.fit = FALSE}, a vector of predictions. If \code{se.fit =
  TRUE}, a list with components

  \item{fit}{ Predictions.}
  \item{se.fit}{ Estimated standard errors.}
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort1Dsmooth}} for computing
  \code{Mort1Dsmooth.object}, \code{\link{predict}}. } 

\examples{# selected data
years <- seq(1850, 2000, 10)
death <- selectHMDdata("Denmark", "Deaths", "Females", 60, years)
exposure <- selectHMDdata("Denmark", "Exposures", "Females", 60, years)
# fit
fit <- Mort1Dsmooth(x=years, y=death, offset=log(exposure),
                    ndx=10, method=3, lambda=10)
# fitted linear predictor at observed times
predict(fit) 
# predict log-rates for each year
plot(fit, "logrates")
newyears <- 1850:2000
pre <- predict(fit, newdata=newyears, type="link", se.fit=TRUE)
lines(newyears, pre$fit, col=4, lty=2, lwd=2)
lines(newyears, pre$fit + 2*pre$se, col=4, lty=3, lwd=2)
lines(newyears, pre$fit - 2*pre$se, col=4, lty=3, lwd=2)
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
