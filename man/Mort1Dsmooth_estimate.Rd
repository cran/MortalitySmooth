\name{Mort1Dsmooth_estimate}

\alias{Mort1Dsmooth_estimate}

\title{ Estimate 1D P-splines for a given lambda}

\description{
  This is an internal function of package \code{MortalitySmooth} which
  estimates coefficients and computes diagnostics for penalized
  B-splines for a given smoothing parameter within the function
  \code{Mort1Dsmooth}.  
}

\usage{
Mort1Dsmooth_estimate(x, y, offset, wei, psi2,
                      B, lambda, DtD, a.init,
                      MON, TOL1, MAX.IT)
}

\arguments{
  \item{x}{ Vector for the abscissa of data. }
  \item{y}{ Vector of counts response. }
  \item{offset}{ Vector with an a priori known component (optional). }
  \item{wei}{ An optional vector of weights to be used in the fitting
              process. } 
  \item{psi2}{ Overdispersion parameter. } 
  \item{B}{ B-splines basis. } 
  \item{lambda}{ Smoothing parameter. } 
  \item{DtD}{ Inner product of the difference matrix. } 
  \item{a.init}{ Vector with the initial coefficients. } 
  \item{MON}{ Logical switch indicating if monitoring is required. } 
  \item{TOL1}{ The tolerance level in the IWLS algorithm. } 
  \item{MAX.IT}{ The maximum number of iterations. }
}
\details{
  Internal function used in \code{\link{Mort1Dsmooth}} for estimating
  coefficients and computing diagnostics. 
}
\value{
  A list with components: 
  \item{a }{ fitted coefficients. }
  \item{h }{ diagonal of the hat-matrix. }
  \item{df}{ effective dimension of used degree of freedom. } 
  \item{aic}{ Akaike's Information Criterion. }
  \item{bic}{ Bayesian Information Criterion. }
  \item{dev}{ Poisson deviance. }
  \item{tol}{ tolerance level. }
  \item{BtWB}{ inner product of basis and weights. }
  \item{P}{ penalty matrix. }
}


\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort1Dsmooth_update}}, \code{\link{Mort1Dsmooth}}. }  

\keyword{smooth}
