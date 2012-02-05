\name{MortSmooth_BWB}
\alias{MortSmooth_BWB}

\title{ Calculates the inner product of a matrix and a sparse weight
  matrix } 

\description{
  This is an internal function of package \code{MortalitySmooth} which
  calculates the inner product of a matrix (from a Kronecker product)
  and a sparse weight matrix within the function \code{Mort2Dsmooth}. It
  construct the LHD of the IWLS algorithm following the idea proposed in
  the Generalized Linear Array Models (see references).
}

\usage{
MortSmooth_BWB(RTBx, RTBy, nbx, nby, W)
}

\arguments{
  \item{RTBx}{ tensors product of B-splines basis for the x-axis. } 
  \item{RTBy}{ tensors product of B-splines basis for the y-axis. } 
  \item{nbx}{ number of B-splines for the x-axis. }
  \item{nby}{ number of B-splines for the y-axis. }
  \item{W}{ matrix of weights. }
}

\details{
  The function employs an arithmetic of arrays which allows to define
  the LHS of the Poisson system of equations as a sequence of nested
  matrix operations. Such way of operating with arrays leads to low
  storage and high speed computation when data are structure as array. 
}

\value{
A matrix of inner product of a matrix and a sparse weight matrix. 
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

\seealso{ \code{\link{Mort2Dsmooth}},
  \code{\link{Mort2Dsmooth_estimate}}, 
  \code{\link{Mort2Dsmooth_update}}. } 

\keyword{ smooth }
\keyword{ regression }
