\name{MortSmooth_BcoefB}
\alias{MortSmooth_BcoefB}

\title{ Multiply a matrix (from a kronecker product) by a matrix }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  compute the linear predictor and the RHS of the IWLS algorithm for 2D
  penalized B-splines within the function \code{Mort2Dsmooth}. 
}

\usage{
MortSmooth_BcoefB(X1, X2, mat)
}

\arguments{
  \item{X1}{ first marginal matrix. }
  \item{X2}{ second marginal matrix. }
  \item{mat}{ matrix to be multiplied. }
}

\details{
Internal function used for computing elements within the
\code{\link{Mort2Dsmooth}} function. 
}

\value{
A matrix that multiplies \code{X1} by \code{mat} and the transpose of
\code{X2}.  
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}}, \code{\link{Mort2Dsmooth_estimate}}. }

\keyword{ smooth }
