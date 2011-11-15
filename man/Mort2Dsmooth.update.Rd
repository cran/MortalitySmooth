\name{Mort2Dsmooth.update}
\alias{Mort2Dsmooth.update}
\title{ Update a 2D Penalized-Poisson Iteration }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  update coefficients for penalized B-splines for two given smoothing
  parameters within the function \code{Mort2Dsmooth}. 
}

\usage{
Mort2Dsmooth.update(x, y, Z, offset, psi2, wei,
                    Bx, By, nbx, nby, RTBx, RTBy, P, a)
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
  \item{P}{ penalty factor. }
  \item{a}{ coefficients (in a matrix). }
}

\details{
  Internal function used for updating a matrix coefficients of the
  B-splines in 2D Poisson penalized B-splines model. 
}

\value{
  A matrix of updated coefficients.
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}}. }

\keyword{smooth}
