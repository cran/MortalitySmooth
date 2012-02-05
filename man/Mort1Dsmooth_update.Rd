\name{Mort1Dsmooth_update}
\alias{Mort1Dsmooth_update}
\title{ Update a 1D Penalized-Poisson Iteration }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  update coefficients for penalized B-splines for a given smoothing
  parameter within the function \code{Mort1Dsmooth}. 
}

\usage{
Mort1Dsmooth_update(x, y, offset, wei, psi2, B, lambdaP, a) 
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
  \item{lambdaP}{ penalty factor, included the smoothing parameter. } 
  \item{a}{ old coefficients. }
}

\details{
  Internal function used for updating coefficients of the B-splines in
  1D Poisson penalized B-splines model. 
}

\value{
  A vector of updated coefficients
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort1Dsmooth}}. }

\keyword{smooth}
