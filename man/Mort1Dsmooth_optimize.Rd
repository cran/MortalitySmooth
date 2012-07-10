\name{Mort1Dsmooth_optimize}
\alias{Mort1Dsmooth_optimize}

\title{ Optimize a 1D Penalized-Poisson IWLS over smoothing parameters }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  optimizes the smoothing parameter for penalized B-splines within the
  function \code{Mort1Dsmooth}.
}

\usage{
Mort1Dsmooth_optimize(x, y, offset, wei,
                      psi2, B, DtD, a.init,
                      MON, TOL1, TOL2, RANGE, MAX.IT,
                      MET)
}

\arguments{
  \item{x}{ vector for the abscissa of data. }
  \item{y}{ vector of counts response. }
  \item{offset}{ vector with an a priori known component (optional). }
  \item{wei}{ an optional vector of weights to be used in the fitting
  process. } 
  \item{psi2}{ an overdispersion parameter used the quasi-likelihood
  approach. } 
  \item{B}{ B-splines basis. }
  \item{DtD}{ inner product of the difference matrix}
  \item{a.init}{ initial coefficients } 
  \item{MON}{ Logical switch indicating if monitoring is required. } 
  \item{TOL1}{ The tolerance level in the IWLS algorithm. }
  \item{TOL2}{ difference between two adjacent smoothing parameters in
    the (pseudo) grid search, log-scale. }
  \item{RANGE}{ range in which smoothing parameter should be searched. }
  \item{MAX.IT}{ the maximum number of iterations }
  \item{MET}{ the method for controlling the amount of smoothing }
}
\details{
  The function aims to find the optimal smoothing parameter within the
  given \code{RANGE} in \code{\link{Mort1Dsmooth}} with \code{method}
  equal to 1 or 2 (BIC and AIC). It employs the function
  \code{cleversearch} from package \code{svcm} in two separate
  steps. First it searches using a rough grid (4 times \code{TOL2}) and
  the median of \code{RANGE} as starting lambda. Afterwards it searches
  in the restricted range around the sub-optimal smoothing parameter,
  using a finer grid defined by \code{TOL2}.

  This procedure allows to find a precise smoothing parameter in an
  efficient way: we do not explore the full range of possible lambda
  values, moving at most one grid step up or down. Furthermore the two
  steps routine reduces the risk of finding sub-optimal smoothing
  parameter.  
}

\author{ Carlo G Camarda }

\references{
Camarda, C. G. (2012). {MortalitySmooth}: An {R} Package for Smoothing
{P}oisson Counts with {P}-Splines. \emph{Journal of Statistical
  Software}. 50, 1-24. \url{http://www.jstatsoft.org/v50/i01/}.
}

\seealso{ \code{\link{Mort1Dsmooth_update}},
  \code{\link{Mort1Dsmooth_estimate}},
  \code{\link{Mort1Dsmooth}}.
}  

\keyword{smooth}
