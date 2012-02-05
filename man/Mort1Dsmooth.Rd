\name{Mort1Dsmooth}

\alias{Mort1Dsmooth}
\alias{print.Mort1Dsmooth}

\title{ Fit One-dimensional Poisson P-splines }

\description{
  Returns an object of class \code{Mort1Dsmooth} which is a P-splines
  smooth of the input data of degree and order fixed by the
  user. Specifically tailored to mortality data.  
}

\usage{
Mort1Dsmooth(x, y, offset, w,
             overdispersion=FALSE,
             ndx = floor(length(x)/5), deg = 3, pord = 2,
             lambda = NULL, df = NULL, method = 1,
             coefstart = NULL,
             control = list())
}

\arguments{
  \item{x}{ a vector with the values of the predictor variable. These must be at least
  \code{2 ndx + 1} of them. } 
  \item{y}{ a vector with a set of counts response variable values. \code{y} must be a
  vector of the same length as \code{x}. }
  \item{offset}{ this can be used to specify an a priori known component
  to be included in the linear predictor during fitting. This should be
  \code{NULL} or a numeric vector of length either one or equal to the
  number of cases. }  
  \item{w}{ an optional vector of weights to be used in the fitting
  process. This should be \code{NULL} or a numeric vector of length
  equal to the number of cases. }
  \item{overdispersion}{ logical on the accounting for overdisperion in
  the smoothing parameter selection criterion. See \code{Details}. Default: FALSE. } 
  \item{ndx}{ number of internal knots -1. Default: \code{floor(length(x)/5)}. } 
  \item{deg}{ degree of the B-splines. Default: 3. } 
  \item{pord}{ order of differences. Default: 2. } 
  \item{lambda}{ smoothing parameter (optional). } 
  \item{df}{ a number which specifies the degrees of freedom (optional). } 
  \item{method}{ the method for controlling the amount of
  smoothing. \code{method = 1} (default) adjusts the smoothing parameter
  so that the BIC is minimized. \code{method = 2} adjusts \code{lambda}
  so that the AIC is minimized. \code{method = 3} uses the value
  supplied for \code{lambda}. \code{method = 4} adjusts \code{lambda} so
  that the degrees of freedom is equal to the supplied \code{df}.}
  \item{coefstart}{ an optional vector of starting coefficients. }
  \item{control}{ a list of control parameters. See \code{Details}. } 
}


\details{
  The method fits a P-spline model with equally-spaced B-splines along
  \code{x}. The response variables must be Poisson distributed counts,
  though overdisperion can be accounted. Offset can be provided,
  otherwise the default is that all weights are one.

  The function is specifically tailored to smooth mortality data in
  one-dimensional setting. In such case the argument \code{x} would be
  either the ages or the years under study. Death counts will be the
  argument \code{y}. In a Poisson regression setting applied to actual
  death counts the \code{offset} will be the logarithm of the exposure
  population. See example below.

  The function can obviously account for zero counts and definite
  offset. In a mortality context, the user can apply the function to
  data with zero deaths, but it has to take care that no exposures are equal
  to zero, i.e. offset equal to minus infinitive. In this last case, the
  argument \code{w} can help. The user would need to set weights equal
  to zero when exposures are equal to zero leading to interpolation of
  the data. See example below.

  Regardless the presence of exposures equal to zero, the argument
  \code{w} can also be used for extrapolation and interpolation of the
  data. Nevertheless see the function
  \code{\link{predict.Mort1Dsmooth}} for a more comprehensive way to
  forecast mortality rates.
  
  The method produces results similar to function \code{smooth.spline},
  but the smoothing function is a B-spline smooth with discrete
  penalization directly on the differences of the B-splines
  coefficients. The user can set the order of difference, the degree of
  the B-splines and number of them. Nevertheless, the smoothing
  parameter \code{lambda} is mainly used to tune the smoothness/model
  fidelity of the fitted values.

  The range in which lambda is searched is given in \code{control} -
  \code{RANGE}. Though it can be modified, the default values are
  suitable for most of the application.  
  
  There are \code{\link{print.Mort1Dsmooth}},
  \code{\link{summary.Mort1Dsmooth}}, \code{\link{plot.Mort1Dsmooth}},
  \code{\link{predict.Mort1Dsmooth}} and
  \code{\link{residuals.Mort1Dsmooth}} methods available for this
  function.  

  Four methods for optimizing the smoothing parameter are available. The BIC
  is set as default. Minimization of the AIC is also possible. BIC will
  give always smoother outcomes with respect to AIC, especially for
  large sample size. Alternatively the user can directly provide the 
  smoothing parameter (\code{method=3}) or the degree of freedom to be
  used in the model (\code{method=4}). Note that \code{Mort1Dsmooth} 
  uses approximated degree of freedom, therefore \code{method=4} will
  produce fitted values with degree of freedom only similar to the one
  provided in \code{df}. The tolerance level can be set via
  \code{control} - \code{TOL2}. 

  Note that the 'ultimate' smoothing with very large lambda will
  approach to a polynomial of degree \code{pord}. 

  Starting coeffients for the B-spline basis can be provided by the
  user. This feature can be useful when a grid-search is manually
  performed by the user. 
  
  The argument \code{overdispersion} can be set to \code{TRUE} when possible presence of
  over(under)dispersion needs to be considered in the selection of the
  smoothing parameter. Mortality data often present overdispersion 
  also known, in demography, as heterogeneity. Duplicates in insurance
  data can lead to overdispersed data, too. Smoothing parameter
  selection may be affected by this phenomenon. When
  \code{overdispersion=TRUE}, the function uses a penalized
  quasi-likelihood method for including an overdisperion parameter
  (\code{psi2}) in the fitting procedure. With this approach expected
  values are assumed equal to the variance multiplied by the parameter
  \code{psi2}. See reference. Note that the inclusion of the
  overdisperion parameter within the estimation might lead to select
  higher lambda, leading to smoother outcomes. When
  \code{overdispersion=FALSE} (default value) or \code{method=3} or
  \code{method=4}, \code{psi2} is estimated after the smoothing
  parameter have been employed. Overdispersion parameter considerably larger
  (smaller) than 1 may be a sign of overdispersion (underdispersion). 
  
  The \code{control} argument is a list that can supply any of the
  following components: 

  \code{MON}: Logical. If \code{TRUE} tracing information on the
  progress of the fitting is produced. Default: \code{FALSE}. 

  \code{TOL1}: The absolute convergence tolerance for each completed
  scoring algorithm. Default: 1e-06. 

  \code{TOL2}: Difference between two adjacent smoothing parameters in
  the (pseudo) grid search, log-scale. Useful only when \code{method} is
  equal to 1, 2 or 4. Default: 0.5.
  
  \code{RANGE}: Range of smoothing parameters in which the grid-search
  is applied, commonly taken in log-scale. Default: [10^-4 ; 10^6].

  \code{MAX.IT}: The maximum number of iterations for each completed
  scoring algorithm. Default: 50.
  
  The arguments \code{MON}, \code{TOL1} and \code{MAX.IT} are kept
  during all the grid search when \code{method} is equal to 1,
  2 or 4. Function \code{cleversearch} from package
  \code{svcm} is employed to speed the grid search. See
  \code{\link{Mort1Dsmooth_optimize}} for details.
 
}


\value{
  An object of the class \code{Mort1Dsmooth} with components:

  \item{coefficients }{ vector of fitted (penalized) B-splines coefficients. } 
  \item{residuals }{ the deviance residuals. }
  \item{fitted.values }{ vector of fitted counts. }
  \item{linear.predictor }{ vector of fitted linear predictor. }
  \item{logmortality }{ fitted mortality rates in log-scale. }
  \item{lev }{ diagonal of the hat-matrix. }
  \item{df }{ effective dimension. }
  \item{deviance }{ Poisson Deviance. }
  \item{aic }{ Akaike's Information Criterion. }
  \item{bic }{ Bayesian Information Criterion. }
  \item{psi2}{ overdispersion parameter. }
  \item{lambda}{ the selected (given) smoothing parameter lambda. } 
  \item{call }{ the matched call. }
  \item{n }{ number of observations. }
  \item{tolerance }{ the used tolerance level. }
  \item{B }{ the B-splines basis. }
  \item{ndx }{ the number of internal knots -1. }
  \item{deg }{ degree of the B-splines. }
  \item{pord }{ order of difference. }
  \item{x }{ values of the predictor variable. }
  \item{y}{ set of counts response variable values. }
  \item{offset}{ vector of the offset. }
  \item{w }{ vector of weights used in the model. }
}

\references{ Eilers and Marx (1996). Flexible Smoothing with B-splines
  and Penalties. Statistical Science. Vol. 11. 89-121. } 

\author{ Carlo G Camarda }

\seealso{ \code{\link{predict.Mort1Dsmooth}},
  \code{\link{plot.Mort1Dsmooth}}, \code{\link{Mort1Dsmooth_optimize}}. }

\examples{
## selected data
years <- 1950:2006
death <- selectHMDdata("Japan", "Deaths", "Females",
                       ages = 80, years = years)
exposure <- selectHMDdata("Japan", "Exposures", "Females",
                          ages = 80, years = years)
## various fits
## default using Bayesian Information Criterion
fitBIC <- Mort1Dsmooth(x=years, y=death,
                       offset=log(exposure))
fitBIC
summary(fitBIC)
## subjective choice of the smoothing parameter lambda
fitLAM <- Mort1Dsmooth(x=years, y=death,
                       offset=log(exposure),
                       method=3, lambda=10000)
## plot
plot(years, log(death/exposure),
main="Mortality rates, log-scale.
      Japanese females, age 80, 1950:2006")
lines(years, fitBIC$logmortality, col=2, lwd=2)
lines(years, fitLAM$logmortality, col=3, lwd=2)
legend("topright", c("Actual", "BIC", "lambda=10000"),
       col=1:3, lwd=c(1,2,2), lty=c(-1,1,1),
       pch=c(1,-1,-1))

## about Extra-Poisson variation (overdispersion)
## checking the presence of overdispersion
fitBIC$psi2 # quite larger than 1
## fitting accounting for overdispersion
fitBICover <- Mort1Dsmooth(x=years, y=death,
                           offset=log(exposure),
                           overdispersion=TRUE)
## difference in the selected smoothing parameters
fitBIC$lambda;fitBICover$lambda
## plotting both situations
plot(fitBICover)
lines(years, fitBIC$logmortality, col=4, lwd=2, lty=2) 

## interpolation, better predict.Mort1Dsmooth
exposure.int <- exposure
exposure.int[20:40] <- NA
## Wrong example (don't run)
## Mort1Dsmooth(x=years, y=death,
##              offset=log(exposure.int))
## using the argument w, with warning
w.int <- rep(1, length(years))
w.int[20:40] <- 0
fit1Dint <- Mort1Dsmooth(x=years, y=death,
                         offset=log(exposure.int), w=w.int)
plot(fit1Dint)

## extrapolation, better using predict.Mort1Dsmooth
years.ext <- seq(years[1], years[length(years)]+20)
exposure.ext <- c(exposure, rep(NA, 20))
death.ext <- c(death, rep(NA, 20))
## Wrong example (don't run)
## Mort1Dsmooth(x=years.ext, y=death.ext,
##              offset=log(exposure.ext))
## using the argument w, with warning
w.ext <- c(rep(1, length(years)), rep(0, 20))
## using the optimal lambda from fitBICover
fit1Dext <- Mort1Dsmooth(x=years.ext, y=death.ext,
                         offset=log(exposure.ext), w=w.ext,
                         method=3, lambda=fitBICover$lambda)
plot(fit1Dext)
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
