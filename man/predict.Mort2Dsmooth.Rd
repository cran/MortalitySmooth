\name{predict.Mort2Dsmooth}
\alias{predict.Mort2Dsmooth}

\title{ Predict Method for 2D P-splines Fits }

\description{
  Obtains predictions and optionally estimates standard errors of those
  predictions from a fitted \code{Mort2Dsmooth} object. 
}

\usage{
## S3 method for class 'Mort2Dsmooth'
\method{predict}{Mort2Dsmooth}(object, newdata = NULL, 
        type = c("link", "response"), 
        se.fit = FALSE, ...)
}

\arguments{
  \item{object}{ An object of class "Mort2Dsmooth", usually, a result of
  a call to \code{Mort2Dsmooth}. } 
  \item{newdata}{ Optionally, a list in which to look for \code{x}
  and/or \code{y} with which to predict. If omitted, the fitted linear
  predictors are used. } 
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
  
  The user can provide also a single predictor (either \code{x} or
  \code{y}) within the argument \code{newdata}. The name within the list
  \code{newdata} must be named \code{x} and \code{y}.  
  
  Forecast is not possible, therefore \code{newdata} has to be include
  within the range of the original \code{x} and \code{y}. 
}

\value{
  If \code{se.fit = FALSE}, a matrix of predictions. If \code{se.fit =
  TRUE}, a list with components 

  \item{fit }{ A matrix of predictions.}
  \item{se.fit }{ A matrix of estimated standard errors.}
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}} for computing
  \code{Mort2Dsmooth.object}, \code{\link{predict}}. } 

\examples{
# selected data
years <- 1980:2006
ages <- 80:100
death <- selectHMDdata("Denmark", "Deaths", "Females",
                       ages = ages, years = years)
exposure <- selectHMDdata("Denmark", "Exposures", "Females",
                          ages = ages, years = years)
# fit
fit <- Mort2Dsmooth(x=ages, y=years, Z=death, offset=log(exposure),
                    method=3, lambdas=c(100,500))

# predict and computing standard errors
pre <- predict(fit, se.fit=TRUE)

# plotting over ages and years 
# !hard to distinguish between upper and lower confidence bounds
grid. <- expand.grid(x = ages, y = years, gr = 1:2)
grid.$lmx <- c(c(pre$fit - 2*pre$se.fit), c(pre$fit + 2*pre$se.fit))
wireframe(lmx ~ x * y, data = grid., groups = gr,
          scales = list(arrows = FALSE),
          drape = TRUE, colorkey = TRUE)

# plotting age 90
plot(years, log(death[11,] / exposure[11,]),
     main="Mortality rates, log-scale. Danish females, age 90, 1980:2006")
lines(years, pre$fit[11,], lwd=2, col=2)
lines(years, pre$fit[11,] + 2*pre$se.fit[11,], lwd=2, col=2, lty=2)
lines(years, pre$fit[11,] - 2*pre$se.fit[11,], lwd=2, col=2, lty=2)

# compute log-death rates for each calendar month and calendar ages
newyears12 <- seq(1990, 2000, length=11*11)
newages12 <- seq(90, 100, length=11*11)
newdata12 <- list(x=newages12, y=newyears12)
pre12 <- predict(fit, newdata=newdata12, se.fit=TRUE)

# death rates in June 1995 at age 95.5
which.age <- which(newages12==95.5)
which.year <- which(newyears12==1995.5)
exp(pre12$fit[which.age, which.year] +
    c(-2*pre12$se.fit[which.age, which.year], 0, 2*pre12$se.fit[which.age, which.year]))
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
