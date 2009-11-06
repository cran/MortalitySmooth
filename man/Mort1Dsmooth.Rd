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
             control = list())
}

\arguments{
  \item{x}{ Values of the predictor variable. These must be at least
  \code{2 ndx + 1} of them. } 
  \item{y}{ Set of counts response variable values. \code{y} must be a
  vector of the same length as \code{x}}  
  \item{offset}{ This can be used to specify an a priori known component
  to be included in the linear predictor during fitting. This should be
  \code{NULL} or a numeric vector of length either one or equal to the
  number of cases. }  
  \item{w}{ An optional vector of weights to be used in the fitting
  process. This should be \code{NULL} or a numeric vector of length
  equal to the number of cases. }
  \item{overdispersion}{ Logical on the accounting for overdisperion in
  the smoothing parameter selection criterion. See
  \code{Details}. Default: FALSE. } 
  \item{ndx}{ Number of internal knots -1. Default:
  \code{floor(length(x)/5)}. } 
  \item{deg}{ Degree of the B-splines. Default: 3. } 
  \item{pord}{ Order of differences. Default: 2. } 
  \item{lambda}{ Smoothing parameter (optional). } 
  \item{df}{ A number which specifies the degrees of freedom
  (optional). } 
  \item{method}{ The method for controlling the amount of
  smoothing. \code{method = 1} (default) adjusts the smoothing parameter
  so that the BIC is minimized. \code{method = 2} adjusts \code{lambda}
  so that the AIC is minimized. \code{method = 3} uses the value
  supplied for \code{lambda}. \code{method = 4} adjusts \code{lambda} so
  that the degrees of freedom is equal to the supplied \code{df}.} 
  \item{control}{ A list of control parameters. See \code{Details}. } 
}


\details{
  The method fits a P-spline model with equally-spaced B-splines along
  \code{x}. The response variables must be Poisson distributed counts,
  though overdisperion can be accounted. Offset can be provided,
  otherwise the default is that all weights are one.

  The method produces results similar to function \code{smooth.spline},
  but the smoothing function is a B-spline smooth with discrete
  penalization directly on the differences of the B-splines
  coefficients. The user can set the order of difference, the degree of
  the B-splines and number of them. Nevertheless, the smoothing
  parameter \code{lambda} is mainly used to tune the smoothness/model
  fidelity of the fitted values.

  There are \code{\link{print.Mort1Dsmooth}},
  \code{\link{summary.Mort1Dsmooth}}, \code{\link{plot.Mort1Dsmooth}}
  \code{\link{predict.Mort1Dsmooth}} and
  \code{\link{residuals.Mort1Dsmooth}} methods available for this
  function. 

  Four methods are possible and optimal smoothing parameter based on BIC
  is set as default. Minimization ofthe AIC is also possible. BIC will
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

  The argument \code{overdispersion} can be set to \code{TRUE} when
  smoothing parameter selection has to account for possible presence of
  over(under)dispersion. Mortality data often present overdispersion
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
  parameter have been employed. Overdispersion parameter larger
  (smaller) than 1 may be a sign of overdispersion (underdispersion). 
  
  The \code{control} argument is a list that can supply any of the
  following components: 

  \code{MON}: Logical. If \code{TRUE} tracing information on the
  progress of the fitting is produced. Default: \code{FALSE}. 

  \code{TOL1}: The absolute convergence tolerance. Default: 1e-06. 

  \code{TOL2}: Difference between two adjacent smoothing parameters in
  the (pseudo) grid search, log-scale. Useful only when \code{method} is
  equal to 1, 2 or 4. Default: 0.1. 

  \code{MAX.IT}: The maximum number of iterations. Default: 50.
  
  The arguments \code{MON}, \code{TOL1} and \code{MAX.IT} are kept
  during all the (pseudo) grid search when \code{method} is equal to 1,
  2 or 4. Function \code{cleversearch} from package
  \code{svcm} is employed to speed the grid search.   
    
  The function is specifically tailored to smooth mortality data in
  one-dimensional setting. In such case the argument \code{x} would be
  either the ages or the years under study. Death counts will be the
  argument \code{y}. In a Poisson regression setting applied to actual
  death counts the \code{offset} will be the logarithm of the exposure
  population. See example below. 
}


\value{
  An object of the class \code{Mort1Dsmooth} with components:

  \item{coefficients }{ vector of fitted (penalized) B-splines
  coefficients. } 
  \item{residuals }{ the deviance residuals. }
  \item{fitted.values }{ vector of fitted counts. }
  \item{linear.predictor }{ vector of fitted linear predictor. }
  \item{leverage }{ diagonal of the hat-matrix. }
  \item{df }{ effective dimension. }
  \item{deviance }{ Poisson Deviance. }
  \item{aic }{ Akaike's Information Criterion. }
  \item{bic }{ Bayesian Information Criterion. }
  \item{psi2}{ Overdispersion parameter. }
  \item{lambda}{ the selected (given) smoothing parameter lambda. } 
  \item{call }{ the matched call. }
  \item{n }{ number of observations. }
  \item{tolerance }{ the used tolerance level. }
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
  \code{\link{plot.Mort1Dsmooth}}. }

\examples{
# selected data
years <- 1950:2006
death <- selectHMDdata("Japan", "Deaths", "Females",
                       ages = 80, years = years)
exposure <- selectHMDdata("Japan", "Exposures", "Females",
                          ages = 80, years = years)
# various fits
# default using Bayesian Information Criterion
fitBIC <- Mort1Dsmooth(x=years, y=death, offset=log(exposure))
fitBIC
summary(fitBIC)
# subjective choice of the smoothing parameter lambda
fitLAM <- Mort1Dsmooth(x=years, y=death, offset=log(exposure),
                       lambda=10000, method=3)
# plot
plot(years, log(death/exposure),
main="Mortality rates, log-scale. Japanese females, age 80, 1950:2006")
lines(years, log(fitted(fitBIC)/exposure), col=2, lwd=2)
lines(years, log(fitted(fitLAM)/exposure), col=3, lwd=2)
legend("topright", c("Actual", "BIC", "lambda=10000"),
       col=1:3, lwd=c(1,2,2), lty=c(-1,1,1),
       pch=c(1,-1,-1))

# about Extra-Poisson variation (overdispersion)
# checking the presence of overdispersion
fitBIC$psi2 # quite larger than 1
# fitting accounting for overdispersion
fitBICover <- Mort1Dsmooth(x=years, y=death, offset=log(exposure),
                           overdispersion=TRUE)
# difference in the selected smoothing parameters
fitBIC$lambda;fitBICover$lambda
# plotting both situations
plot(fitBICover)
lines(years, log(fitBIC$fitted) - fitBIC$offset, col=4, lwd=2, lty=2) 
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
