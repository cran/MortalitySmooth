\name{predict.Mort2Dsmooth}
\alias{predict.Mort2Dsmooth}

\title{ Predict Method for 2D P-splines Fits }

\description{
  Obtains predictions, forecasts and optionally estimated standard
  errors of those predictions from a fitted \code{Mort2Dsmooth} object.  
}

\usage{
\method{predict}{Mort2Dsmooth}(object, newdata = NULL, 
        type = c("link", "response"), 
        se.fit = FALSE, ...)
}

\arguments{
  \item{object}{ an object of class "Mort2Dsmooth", usually, a result of
  a call to \code{Mort2Dsmooth}. } 
  \item{newdata}{ optionally, a list in which to look for \code{x}
  and/or \code{y} with which to predict. If omitted, the fitted linear
  predictors are used. } 
  \item{type}{ the type of prediction required. The default ("link") is
  on the scale of the linear predictors; the alternative "response" is
  on the scale of the response variable. } 
  \item{se.fit}{ logical switch indicating if standard errors are
    required. Default: \code{FALSE}. }
  \item{\ldots}{ other predict parameters to passed to
  \code{predict}. Not in used. }
}

\details{
  If \code{newdata} is omitted the predictions are based on the data
  used for the fit. Note that, in common with other prediction
  functions, any offset supplied as an argument is always ignored when
  predicting, unlike offsets specified in modelling. 
  
  The user can provide also a single predictor (either \code{x} or
  \code{y}) within the argument \code{newdata}. The name within the list
  \code{newdata} must be named \code{x} and \code{y}.
  
  Forecast is possible, therefore \code{newdata} can include
  values out of the range of the original \code{x} and \code{y}. See
  example below. 

  Interpolation is also feasible. See example below.
}

\value{
  If \code{se.fit = FALSE}, a matrix of predictions. If \code{se.fit =
  TRUE}, a list with components:

  \item{fit }{ a matrix of predictions.}
  \item{se.fit }{ a matrix of estimated standard errors.}
}

\author{ Carlo G Camarda }

\references{
Camarda, C. G. (2012). {MortalitySmooth}: An {R} Package for Smoothing
{P}oisson Counts with {P}-Splines. \emph{Journal of Statistical
Software}. 50, 1-24. \url{http://www.jstatsoft.org/v50/i01/}.

Currie, I. D., M. Durban, and P. H. C. Eilers (2004). Smoothing and
forecasting mortality rates. \emph{Statistical Modelling}. 4, 279-298. 
}

\seealso{ \code{\link{Mort2Dsmooth}} for computing
  \code{Mort2Dsmooth.object}. } 

\examples{
## computing confidence intervals
## selected data
years <- 1980:2006
ages <- 80:100
death <- selectHMDdata("Denmark", "Deaths", "Females",
                       ages = ages, years = years)
exposure <- selectHMDdata("Denmark", "Exposures", "Females",
                          ages = ages, years = years)
## fit
fit <- Mort2Dsmooth(x=ages, y=years, Z=death,
                    offset=log(exposure),
                    method=3, lambdas=c(100,500))

## predict and computing standard errors
pre <- predict(fit, se.fit=TRUE)

## plotting over ages and years 95% confidence intervals 
## !hard to distinguish between upper
## and lower confidence bounds
grid. <- expand.grid(x = ages, y = years, gr = 1:2)
grid.$lmx <- c(c(pre$fit - 2*pre$se.fit),
               c(pre$fit + 2*pre$se.fit))
wireframe(lmx ~ x * y, data = grid., groups = gr,
          scales = list(arrows = FALSE),
          drape = TRUE, colorkey = TRUE)

## plotting age 80 (first row)
plot(years, log(death[1,] / exposure[1,]),
     main="Mortality rates, log-scale.
           Danish females, age 80, 1980:2006")
lines(years, pre$fit[1,], lwd=2, col=2)
lines(years, pre$fit[1,] + 2*pre$se.fit[1,],
      lwd=2, col=2, lty=2)
lines(years, pre$fit[1,] - 2*pre$se.fit[1,],
      lwd=2, col=2, lty=2)

## forecasting example
newyears <- 1980:2020
newdata <- list(x=ages, y=newyears) 
pre.for <- predict(fit, newdata=newdata, se.fit=TRUE)

## plot fitted+forecast log-rates for all ages over years
matplot(years, t(log(death/exposure)), pch=1, cex=0.5,
        col=rainbow(length(ages)),
        xlim=range(newyears),
        ylim=range(pre.for$fit),
        ylab="log-mortality")
matlines(newyears, t(pre.for$fit), lty=1, lwd=2,
        col=rainbow(length(ages)))

## plot selected ages over years, with 95% confidence intervals
whiA <- c(1,6,11)
matplot(years, t(log(death[whiA,]/exposure[whiA,])),
        pch=1,
        xlim=range(newyears),
        ylim=c(-3.3, -1.5),
        ylab="log-mortality")
matlines(newyears, t(pre.for$fit[whiA,]), lty=1, lwd=2)
matlines(newyears, t(pre.for$fit[whiA,]+
                     2*pre.for$se.fit[whiA,]), lty=2)
matlines(newyears, t(pre.for$fit[whiA,]-
                     2*pre.for$se.fit[whiA,]), lty=2)


## interpolation example
## compute log-death rates for
## each calendar month and calendar ages
newyears12 <- seq(1990, 2000, length=11*11)
newages12 <- seq(90, 100, length=11*11)
newdata12 <- list(x=newages12, y=newyears12)
pre12 <- predict(fit, newdata=newdata12, se.fit=TRUE)

## death rates in June 1995 at age 95.5
which.age <- which(newages12==95.5)
which.year <- which(newyears12==1995.5)
exp(pre12$fit[which.age, which.year] +
    c(-2*pre12$se.fit[which.age, which.year],
      0, 2*pre12$se.fit[which.age, which.year]))
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
