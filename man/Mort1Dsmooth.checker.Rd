\name{Mort1Dsmooth.checker}
\alias{Mort1Dsmooth.checker}

\title{ Check Arguments for Function Mort1Dsmooth }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  checks whether the arguments given in the function \code{Mort1Dsmooth}
  have proper lengths and suitability. 
}

\usage{
Mort1Dsmooth.checker(x, y, offset, w,
                     overdispersion,
                     ndx, deg, pord, lambda, df,
                     method, control) 
}

\arguments{
  \item{x}{ Values of the predictor variable. } 
  \item{y}{ Set of counts response variable values.  } 
  \item{offset}{ A vector with the offset. } 
  \item{w}{ An optional vector of weights to be used in the fitting
            process. } 
  \item{overdispersion}{ Logical on the accounting for overdisperion. }
  \item{ndx}{ Number of internal knots -1.  } 
  \item{deg}{ Degree of the B-splines.  } 
  \item{pord}{ Order of differences.  } 
  \item{lambda}{ Smoothing parameter. } 
  \item{df}{ A number which specifies the degrees of freedom. } 
  \item{method}{ The method for controlling the amount of smoothing. } 
  \item{control}{ A list of control parameters. } 
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
