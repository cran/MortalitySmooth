\name{Mort1Dsmooth.update}
\alias{Mort1Dsmooth.update}
\title{ Update a 1D Penalized-Poisson Iteration }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  update coefficients for penalized B-splines for a given smoothing
  parameter within the function \code{Mort1Dsmooth}. 
}

\usage{
Mort1Dsmooth.update(x, y, offset, wei, psi2, B, lambdaP, a) 
}

\arguments{
  \item{x}{ Vector for the abscissa of data }
  \item{y}{ Vector of counts response }
  \item{offset}{ Vector with an a priori known component (optional) }
  \item{wei}{ An optional vector of weights to be used in the fitting
  process } 
  \item{psi2}{ An overdispersion parameter used the quasi-likelihood
  approach } 
  \item{B}{ B-splines basis }
  \item{lambdaP}{ Penalty factor, included the smoothing parameter } 
  \item{a}{ Old coefficients }
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
