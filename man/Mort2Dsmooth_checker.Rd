\name{Mort2Dsmooth_checker}
\alias{Mort2Dsmooth_checker}

\title{ Check Arguments for Function Mort2Dsmooth }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  checks whether the arguments given in the function \code{Mort2Dsmooth}
  have proper lengths and suitability. 
}

\usage{
Mort2Dsmooth_checker(x, y, Z, offset, W,
                     overdispersion,
                     ndx, deg, pord, lambdas,
                     df, method, coefstart,
                     control) 
}

\arguments{
  \item{x}{ vector for the abscissa of data. }
  \item{y}{ vector for the ordinate of data. }
  \item{Z}{ matrix of counts response. }
  \item{offset}{ matrix with an a priori known component. }
  \item{W}{ an optional matrix of weights to be used in the fitting
  process. } 
  \item{overdispersion}{ logical on the accounting for overdisperion. }
  \item{ndx}{ vector with the number of internal knots -1 for each
  axis. } 
  \item{deg}{ vector with the degree of the B-splines for each axis.  } 
  \item{pord}{ vector with the order of differences for each axis.  }
  \item{lambdas}{ vector with smoothing parameters, possibly one for
  axis. } 
  \item{df}{ a number which specifies the degrees of freedom. }
  \item{method}{ the method for controlling the amount of smoothing. }
  \item{coefstart}{ a matrix with the eventual starting coefficients. }
  \item{control}{ a list of control parameters. }
}

\details{
  Internal function used in \code{\link{Mort2Dsmooth}} for checking its
  arguments. 
}

\value{
  A list with checked elements to be used in
  \code{\link{Mort2Dsmooth}}. 
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}}. }

\keyword{smooth}
