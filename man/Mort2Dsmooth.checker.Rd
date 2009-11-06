\name{Mort2Dsmooth.checker}
\alias{Mort2Dsmooth.checker}

\title{ Check Arguments for Function Mort2Dsmooth }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  checks whether the arguments given in the function \code{Mort2Dsmooth}
  have proper lengths and suitability. 
}

\usage{
Mort2Dsmooth.checker(x, y, Z, offset, W,
                     overdispersion,
                     ndx, deg, pord, lambdas,
                     df, method, control) 
}

\arguments{
  \item{x}{ Vector for the abscissa of data }
  \item{y}{ Vector for the ordinate of data }
  \item{Z}{ Matrix of counts response }
  \item{offset}{ Matrix with an a priori known component }
  \item{W}{ An optional matrix of weights to be used in the fitting
  process } 
  \item{overdispersion}{ Logical on the accounting for overdisperion. }
  \item{ndx}{ Vector with the number of internal knots -1 for each
  axis. } 
  \item{deg}{ Vector with the degree of the B-splines for each axis.  } 
  \item{pord}{ Vector with the order of differences for each axis.  }
  \item{lambdas}{ Vector with smoothing parameters, possibly one for
  axis. } 
  \item{df}{ A number which specifies the degrees of freedom. }
  \item{method}{ The method for controlling the amount of smoothing. } 
  \item{control}{ A list of control parameters. }
}

\details{
  Internal function used in \code{\link{Mort2Dsmooth}} for checking its
  arguments. 
}

\value{
  A list with checked elements to be used in \code{\link{Mort2Dsmooth}} 
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}}. }

\keyword{smooth}
