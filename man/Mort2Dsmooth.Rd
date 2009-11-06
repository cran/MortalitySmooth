\name{Mort2Dsmooth}

\alias{Mort2Dsmooth}
\alias{print.Mort2Dsmooth}

\title{ Fit Two-dimensional Poisson P-splines }

\description{
  Returns an object of class \code{Mort2Dsmooth} which is a
  two-dimensional P-splines smooth of the input data of degree and order
  fixed by the user. Specifically tailored to mortality data. 
}

\usage{
Mort2Dsmooth(x, y, Z, offset, W, overdispersion=FALSE, 
             ndx = c(floor(length(x)/5), floor(length(y)/5)), 
             deg = c(3, 3), pord = c(2, 2), 
             lambdas = NULL, df = NULL, method = 1, 
             control = list())
}

\arguments{
  \item{x}{ Vector for the abscissa of data. These must be at least
  \code{2 ndx[1] + 1} of them. } 
  \item{y}{ Vector for the ordinate of data. These must be at least
  \code{2 ndx[2] + 1} of them. } 
  \item{Z}{ Matrix of counts response. Dimensions of \code{Z} should
  correspond to the length of \code{x} and \code{y}, respectively. } 
  \item{offset}{ Matrix with an a priori known component to be included
  in the linear predictor during fitting. This should be \code{NULL} or
  a numeric matrix of dimensions of \code{Z} or a single numeric
  value. } 
  \item{W}{ An optional matrix of weights to be used in the fitting
  process. This should be \code{NULL} or a numeric matrix of dimensions
  of \code{Z} or a single numeric value. }
  \item{overdispersion}{ Logical on the accounting for overdisperion in
  the smoothing parameters selection criterion. See
  \code{Details}. Default: FALSE. }	   
  \item{ndx}{ Vector with the number of internal knots -1 for each
  axis. Default: [\code{floor(length(x)/5)},
  \code{floor(length(y)/5)}]. } 
  \item{deg}{ Vector with the degree of the B-splines for each
  axis. Default: [3, 3]. } 
  \item{pord}{ Vector with the order of differences for each
  axis. Default: [2, 2]. } 
  \item{lambdas}{ Vector with smoothing parameters, possibly one for
  axis (optional). } 
  \item{df}{ A number which specifies the degrees of freedom
  (optional). } 
  \item{method}{ The method for controlling the amount of
  smoothing. \code{method = 1} (default) adjusts the two smoothing
  parameters so that the BIC is minimized. \code{method = 2} adjusts
  \code{lambdas} so that the AIC is minimized. \code{method = 3} uses
  the values supplied for \code{lambdas}. Isotropic smoothing is allowed
  in this method. \code{method = 4} adjusts \code{lambdas} so that the
  degrees of freedom is equal to the supplied \code{df}. } 
  \item{control}{ A list of control parameters. See \code{Details}. }
}


\details{
  The method fits a two-dimensional P-spline model with equally-spaced
  B-splines along \code{x} and \code{y}. The response variables must be
  a matrix of Poisson distributed counts. Offset can be provided,
  otherwise the default is that all weights are one. 

  The method produces results from a smoothing function which is the
  Kronecker product of B-spline basis over the two axes and include a
  discrete penalization directly on the differences of the B-splines
  coefficients. The user can set the order of difference, the degree of
  the B-splines and number of them for each of the axis. Nevertheless,
  the smoothing parameters \code{lambdas} are mainly used to tune the 
  smoothness/model fidelity of the fitted values. 
  
  There are \code{\link{print.Mort2Dsmooth}},
  \code{\link{summary.Mort2Dsmooth}}, \code{\link{plot.Mort2Dsmooth}}
  \code{\link{predict.Mort2Dsmooth}} and
  \code{\link{residuals.Mort2Dsmooth}} methods available for this 
  function. 

  Four methods are possible and the minimization of the BIC is set as
  default method for optimizing the two smoothing
  parameters. Minimization of the AIC is also possible. BIC will give
  always smoother outcomes with respect to AIC, especially for large
  sample size. Alternatively the user can directly provide the smoothing
  parameters (\code{method=3}) or the degree of freedom to be used in
  the model (\code{method=4}). If the user provides only a single value
  for the argument \code{lambdas}, isotropic smoothing is applied (with
  warning). Note that \code{Mort2Dsmooth} uses approximated degree of
  freedom, therefore \code{method=4} will produce fitted values with
  degree of freedom only similar to the one provided in \code{df}. The
  tolerance level can be set via \code{control} - \code{TOL2}. 

  Note that the two-dimensional 'ultimate' smoothing with very large
  lambda will approach to a surface which is a product of two polynomial
  of degree \code{pord[1]} and \code{pord[2]}, respectively. In
  particular, when \code{pord=c(2,2)} the 'ultimate'  smoothing is a
  bi-linear surface over \code{x} and \code{y}. 

  The argument \code{overdispersion} can be set to \code{TRUE} when
  smoothing parameters selection has to account for possible presence of 
  over(under)dispersion. Mortality data often present overdispersion
  also known, in demography, as heterogeneity. Duplicates in insurance
  data can lead to overdispersed data, too. Smoothing parameters
  selection may be affected by this phenomenon. When
  \code{overdispersion=TRUE}, the function uses a penalized
  quasi-likelihood method for including an overdisperion parameter
  (\code{psi2}) in the fitting procedure. With this approach expected
  values are assumed equal to the variance multiplied by the parameter
  \code{psi2}. See reference. Note that with overdispersed data both BIC
  and AIC might select higher lambdas, leading to smoother
  outcomes. When \code{overdispersion=FALSE} (default value) or
  \code{method=3} or \code{method=4}, \code{psi2} is estimated after the
  smoothing parameters have been employed. Overdispersion parameter
  larger (smaller) than 1 may be a sign of overdispersion
  (underdispersion). 

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
  either the ages and the argument \code{y} the years under study. The
  matrix of death counts will be the argument \code{Z}. In a Poisson
  regression setting applied to actual death counts the \code{offset}
  will be the logarithm of the matrix of exposure population. See
  example below. 
  
  The inner functions work using an arithmetic of arrays defined as
  Generalized Linear Array Model (GLAM) (see reference). In order to
  avoid construction of large Kronecker product basis from the large
  number of B-splines along the axes, the function profits of the
  special structure of both the data as rectangular array and the model
  matrix as tensor product. It uses sequence of  nested matrix
  operations and this leads to low storage and high speed computation
  within the IWLS algorithm. Moreover, the function do not vectorize the
  whole system keeping the actual two-dimensional array structure within
  the scoring algorithm. 
}

\value{
  An object of the class \code{Mort2Dsmooth} with components:

  \item{coefficients }{ matrix of fitted (penalized) B-splines
  coefficients. } 
  \item{residuals }{ matrix of deviance residuals. }
  \item{fitted.values }{ matrix of fitted counts. }
  \item{linear.predictor }{ matrix of fitted linear predictor. }
  \item{df  }{ effective dimension. }
  \item{dev }{ Poisson Deviance. }
  \item{aic }{ Akaike's Information Criterion. }
  \item{bic }{ Bayesian Information Criterion. }
  \item{psi2 }{ Overdispersion parameter. }
  \item{lambdas}{ the selected (given) smoothing parameters. }
  \item{call }{ the matched call. }
  \item{m }{ length of argument \code{x}. }
  \item{n }{ length of argument \code{y}. }
  \item{tolerance }{ the used tolerance level. }
  \item{lev }{ diagonal of the hat-matrix. }
  \item{ndx }{ the number of internal knots -1 for each axis. }
  \item{deg }{ degree of the B-splines for each axis. }
  \item{pord }{ order of difference for each axis. }
  \item{x}{ vector for the abscissa of data. }
  \item{y}{ vector for the ordinate of data. }
  \item{Z}{ matrix of counts response. }
  \item{offset}{ matrix with an a priori known component. }
  \item{w}{ matrix of weights. }
}

\references{ 
Eilers, P. H. C., I. D. Currie, and M. Durban (2006). Fast and compact
smoothing on large multidimensional grids. \emph{Computational
  Statistics \& Data Analysis} 50, 61-76. 

Currie, I. D., M. Durban, and P. H. C. Eilers (2006). Generalized linear
array models with applications to multidimentional
smoothing. \emph{Journal of the Royal Statistical Society. Series
  B}. 68, 259-280. 
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{predict.Mort2Dsmooth}},
  \code{\link{plot.Mort2Dsmooth}}. } 

\examples{
# selected data
ages <- 50:100
years <- 1950:2006
death <- selectHMDdata("Japan", "Deaths", "Females",
                       ages = ages, years = years) 
exposure <- selectHMDdata("Japan", "Exposures", "Females",
                          ages = ages, years = years)

# fit with BIC
fitBIC <- Mort2Dsmooth(x=ages, y=years, Z=death, offset=log(exposure)) 
fitBIC
summary(fitBIC)

# plot age 50 log death rates (1st row)
plot(years, log(death[1,]/exposure[1,]),
     main="Mortality rates, log-scale. Japanese females, age 50, 1950:2006")
lines(years, log(fitted(fitBIC)[1,]/exposure[1,]), col=2, lwd=2)

# plot over age and years fitted log death rates from fitBIC
grid. <- expand.grid(list(ages=ages, years=years))
grid.$lmx <- c(log(fitted(fitBIC)/exposure))
levelplot(lmx ~ years * ages , grid., 
          at=quantile(grid.$lmx, seq(0,1,length=10)),
          col.regions=rainbow(9))

# about Extra-Poisson variation (overdispersion)
# checking the presence of overdispersion
fitBIC$psi2
# fitting accounting for overdispersion
fitBICover <- Mort2Dsmooth(x=ages, y=years, Z=death,
                           offset=log(exposure),
                           overdispersion=TRUE)
# difference in the selected smoothing parameters
fitBIC$lambdas;fitBICover$lambdas
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
