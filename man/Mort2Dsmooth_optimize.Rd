\name{Mort2Dsmooth_optimize}
\alias{Mort2Dsmooth_optimize}

\title{ Optimize a 2D Penalized-Poisson IWLS over smoothing parameters }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  optimizes the smoothing parameter for penalized B-splines within the
  function \code{Mort2Dsmooth}.
}

\usage{
Mort2Dsmooth_optimize(x, y, Z, offset, wei,
                      psi2, Bx, By, nbx, nby,
                      RTBx, RTBy, Px, Py,
                      a.init,
                      MON, TOL1, TOL2,
                      RANGEx, RANGEy,
                      MAX.IT, MET)
}

\arguments{
  \item{x}{ vector for the abscissa of data. }
  \item{y}{ vector for the ordinate of data. }
  \item{Z}{ matrix of counts response. }
  \item{offset}{ matrix with an a priori known component (optional). }
  \item{wei}{ an optional matrix of weights to be used in the fitting
  process. } 
  \item{psi2}{ overdispersion parameter. }
  \item{Bx}{ B-splines basis for the x-axis. }
  \item{By}{ B-splines basis for the y-axis. }
  \item{nbx}{ number of B-splines for the x-axis. }
  \item{nby}{ number of B-splines for the y-axis. }
  \item{RTBx}{ tensors product of B-splines basis for the x-axis. }
  \item{RTBy}{ tensors product of B-splines basis for the y-axis. }
  \item{Px}{ penalty factor for the x-axis. }
  \item{Py}{ penalty factor for the y-axis. }
  \item{a.init}{ matrix with the initial coefficients. }
  \item{MON}{ Logical switch indicating if monitoring is required. } 
  \item{TOL1}{ The tolerance level in the IWLS algorithm. }
  \item{TOL2}{ difference between two adjacent smoothing parameters in
    the grid search, log-scale. }
  \item{RANGEx}{ range in which smoothing parameter for x should be
    searched. }
  \item{RANGEy}{ range in which smoothing parameter for y should be
  searched. }
  \item{MAX.IT}{ the maximum number of iterations }
  \item{MET}{ the method for controlling the amount of smoothing }
}
\details{
  The function aims to find the optimal smoothing parameters within the
  given \code{RANGEx} and \code{RANGEy} in \code{\link{Mort2Dsmooth}}
  with \code{method} equal to 1 or 2 (BIC and AIC). It employs the
  function \code{cleversearch} from package \code{svcm} in two
  separate steps. First it searches using a rough grid (4 times
  \code{TOL2}) and the median of \code{RANGEx} and \code{RANGEy} as
  starting lambdas. Afterwards it searches in the restricted areas
  around the sub-optimal smoothing parameters, using a finer grid
  defined by \code{TOL2}. 

  This procedure allows to find precise smoothing parameters in an
  efficient way: we do not explore the full ranges of possible lambda
  values, but we optimize each parameter in turn, moving at most one
  grid step up or down. Furthermore the two steps routine reduces the
  risk of finding sub-optimal smoothing parameters.   
  }

\author{ Carlo G Camarda }

\references{
Camarda, C. G. (2012). {MortalitySmooth}: An {R} Package for Smoothing
{P}oisson Counts with {P}-Splines. \emph{Journal of Statistical
Software}. 50, 1-24. \url{http://www.jstatsoft.org/v50/i01/}.
}


\seealso{ \code{\link{Mort2Dsmooth_update}},
  \code{\link{Mort2Dsmooth_estimate}},
  \code{\link{Mort2Dsmooth}}.
}  
\keyword{smooth}
