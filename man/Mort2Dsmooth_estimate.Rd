\name{Mort2Dsmooth_estimate}

\alias{Mort2Dsmooth_estimate}

\title{ Estimate 2D P-splines for two given lambdas}

\description{
  This is an internal function of package \code{MortalitySmooth} which
  estimates coefficients and computes diagnostics for two-dimensional
  penalized B-splines for two given smoothing parameters within the
  function \code{Mort2Dsmooth}. 
}

\usage{
Mort2Dsmooth_estimate(x, y, Z, offset, psi2, wei,
                      Bx, By, nbx, nby, RTBx, RTBy,
                      lambdas, Px, Py, a.init, 
                      MON, TOL1, MAX.IT)
}

\arguments{
  \item{x}{ vector for the abscissa of data. }
  \item{y}{ vector for the ordinate of data. }
  \item{Z}{ matrix of counts response. }
  \item{offset}{ matrix with an a priori known component (optional). } 
  \item{psi2}{ overdispersion parameter. }
  \item{wei}{ an optional matrix of weights to be used in the fitting
  process. } 
  \item{Bx}{ B-splines basis for the x-axis. }
  \item{By}{ B-splines basis for the y-axis. }
  \item{nbx}{ number of B-splines for the x-axis. }
  \item{nby}{ number of B-splines for the y-axis. }
  \item{RTBx}{ tensors product of B-splines basis for the x-axis. }
  \item{RTBy}{ tensors product of B-splines basis for the y-axis. }
  \item{lambdas}{ vector with the two smoothing parameters. }
  \item{Px}{ penalty factor for the x-axis. }
  \item{Py}{ penalty factor for the y-axis. }
  \item{a.init}{ matrix with the initial coefficients. }
  \item{MON}{ logical switch indicating if monitoring is required. }
  \item{TOL1}{ the tolerance level in the IWLS algorithm. }
  \item{MAX.IT}{ the maximum number of iterations. }
}

\details{
  Internal function used in \code{\link{Mort2Dsmooth}} for estimating
  coefficients and computing diagnostics. 
}

\value{
  A list with components: 
  \item{a }{ fitted coefficients (in a matrix). }
  \item{h }{ diagonal of the hat-matrix. }
  \item{df}{ effective dimension of used degree of freedom. } 
  \item{aic}{ Akaike's Information Criterion. }
  \item{bic}{ Bayesian Information Criterion. }
  \item{dev}{ Poisson deviance. }
  \item{tol}{ tolerance level. }
  \item{BWB}{ inner product of basis and weights. }
  \item{P}{ penalty matrix. }
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth_update}},
  \code{\link{Mort2Dsmooth}}. } 

\keyword{smooth}
