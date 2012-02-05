\name{Mort1Dsmooth_checker}
\alias{Mort1Dsmooth_checker}

\title{ Check Arguments for Function Mort1Dsmooth }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  checks whether the arguments given in the function \code{Mort1Dsmooth}
  have proper lengths and suitability. 
}

\usage{
Mort1Dsmooth_checker(x, y, offset,
                     w, overdispersion,
                     ndx, deg, pord, lambda,
                     df, method, coefstart,
                     control)
}

\arguments{
  \item{x}{ values of the predictor variable. } 
  \item{y}{ set of counts response variable values.  } 
  \item{offset}{ a vector with the offset. } 
  \item{w}{ an optional vector of weights to be used in the fitting
            process. } 
  \item{overdispersion}{ logical on the accounting for overdisperion. }
  \item{ndx}{ number of internal knots -1.  } 
  \item{deg}{ degree of the B-splines.  } 
  \item{pord}{ order of differences.  } 
  \item{lambda}{ smoothing parameter. } 
  \item{df}{ a number which specifies the degrees of freedom. } 
  \item{method}{ the method for controlling the amount of smoothing. }
  \item{coefstart}{ a vector with the eventual starting coefficients. }
  \item{control}{ a list of control parameters. } 
}

\details{
  Internal function used in \code{\link{Mort1Dsmooth}} for checking its
  arguments. 
}

\value{
  A list with checked elements to be used in \code{\link{Mort1Dsmooth}} 
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort1Dsmooth}} }

\keyword{smooth}
