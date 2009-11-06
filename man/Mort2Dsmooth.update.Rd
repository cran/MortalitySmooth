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
  \item{x}{ Vector for the abscissa of data. }
  \item{y}{ Vector for the ordinate of data. }
  \item{Z}{ Matrix of counts response. }
  \item{offset}{ Matrix with an a priori known component (optional). }
  \item{wei}{ An optional matrix of weights to be used in the fitting
  process. } 
  \item{psi2}{ Overdispersion parameter. }
  \item{Bx}{ B-splines basis for the x-axis. }
  \item{By}{ B-splines basis for the y-axis. }
  \item{nbx}{ Number of B-splines for the x-axis. }
  \item{nby}{ Number of B-splines for the y-axis. }
  \item{RTBx}{ Tensors product of B-splines basis for the x-axis. }
  \item{RTBy}{ Tensors product of B-splines basis for the y-axis. }
  \item{P}{ Penalty factor. }
  \item{a}{ Coefficients (in a matrix). }
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
