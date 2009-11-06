\name{MortSmooth.bbase}
\alias{MortSmooth.bbase}

\title{ Construct B-spline basis }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  creates equally-spaced B-splines basis over an abscissa of data within
  the function \code{Mort1Dsmooth}. 
}

\usage{
MortSmooth.bbase(x, xl, xr, ndx, deg)
}

\arguments{
  \item{x}{ Vector for the abscissa of data. }
  \item{xl}{ Left boundary. }
  \item{xr}{ Right boundary. }
  \item{ndx}{ Number of internal knots minus one or number of internal
  intervals. } 
  \item{deg}{ Degree of the splines. }
}

\details{
  The function reproduce an algorithm presented by Eilers and Marx
  (2004) using differences of truncated power functions (see
  \code{\link{MortSmooth.tpower}}). The final matrix has a single
  B-spline for each of the [\code{ndx} + \code{deg}] columns. The number
  of rows is equal to the length of \code{x}. 
}

\value{
  A matrix containing equally-spaced B-splines of degree \code{deg}
  along \code{x} for each column.  
}

\references{ Eilers and Marx (2004). Splines, Knots, and
  Penalties. Unpublished manuscript.} 

\author{ Carlo G Camarda }

\seealso{ \code{\link{MortSmooth.tpower}}. }

\examples{
x <- seq(0,1,length=100)
# B-splines basis of degree 1
B1 <- MortSmooth.bbase(x=x, xl=min(x), xr=max(x), ndx=10, deg=1)
matplot(x, B1, t="l", main="B-splines basis of degree 1")
# B-splines basis of degree 3
B3 <- MortSmooth.bbase(x=x, xl=min(x), xr=max(x), ndx=10, deg=3)
matplot(x, B3, t="l", main="B-splines basis of degree 3")
}

\keyword{smooth}
