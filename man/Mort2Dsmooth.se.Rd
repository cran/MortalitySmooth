\name{Mort2Dsmooth.se}
\alias{Mort2Dsmooth.se}

\title{ Compute a 2D standard errors }

\description{
   This is an internal function of package \code{MortalitySmooth} which
   calculates the inner product of a matrix (from a Kronecker product)
   and a sparse weight matrix in order to obtain standard errors. It
   uses the same idea employed in \code{MortSmooth.BWB} and the elements
   after the IWLS converged, including the penalty term. 
}

\usage{
Mort2Dsmooth.se(RTBx, RTBy, nbx, nby, BWB.P1)
}


\arguments{
  \item{RTBx}{ Tensors product of B-splines basis for the x-axis. }
  \item{RTBy}{ Tensors product of B-splines basis for the y-axis. }
  \item{nbx}{ Number of B-splines for the x-axis. }
  \item{nby}{ Number of B-splines for the y-axis. }
  \item{BWB.P1}{ Inverse of the LHS of the Poisson system of equations. } 
}

\details{
  This function is only used within \code{\link{predict.Mort2Dsmooth}}
  when standard errors are required. The arguments \code{BWB.P1} is the
  LHS after convergence is reached and smoothing parameter selected. The
  standard errors as given in the function are computed for the linear
  predictor term and simple computation is needed to obtain standard
  errors for the Poisson counts. Anyway
  \code{\link{predict.Mort2Dsmooth}} takes care of such differences.  
  
  The Generalized Linear Array Models setting is explained in the
  reference in \code{\link{MortSmooth.BWB}} and
  \code{\link{Mort2Dsmooth}}. 
}

\value{
A matrix of standard errors for the linear predictor term.
}

\author{ Carlo G Camarda }


\seealso{ \code{\link{Mort2Dsmooth}}, \code{\link{MortSmooth.BWB}},
  \code{\link{predict.Mort2Dsmooth}}. } 

\keyword{ smooth }
\keyword{ regression }
