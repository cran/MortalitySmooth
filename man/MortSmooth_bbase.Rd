\name{MortSmooth_bbase}
\alias{MortSmooth_bbase}

\title{ Construct B-spline basis }

\description{
  This is an internal function of package \code{MortalitySmooth} which
  creates equally-spaced B-splines basis over an abscissa of data within
  the function \code{Mort1Dsmooth}. 
}

\usage{
MortSmooth_bbase(x, xl, xr, ndx, deg)
}

\arguments{
  \item{x}{ vector for the abscissa of data. }
  \item{xl}{ left boundary. }
  \item{xr}{ right boundary. }
  \item{ndx}{ number of internal knots minus one or number of internal
  intervals. } 
  \item{deg}{ degree of the splines. }
}

\details{
  The function reproduce an algorithm presented by Eilers and Marx
  (2010) using differences of truncated power functions (see
  \code{\link{MortSmooth_tpower}}). The final matrix has a single
  B-spline for each of the [\code{ndx} + \code{deg}] columns. The number
  of rows is equal to the length of \code{x}.

  The function differs from \code{bs} in the package \code{splines}
  since it automatically constructed B-splines with identical
  shape. This would allow a simple interpretation of coefficients and
  application of simple differencing.
}

\value{
  A matrix containing equally-spaced B-splines of degree \code{deg}
  along \code{x} for each column.  
}

\author{ Carlo G Camarda }

\references{
Eilers P. H. C. and B. D. Marx (2010). Splines, Knots, and
Penalties. \emph{Wiley Interdisciplinary Reviews: Computational
  Statistics}. 2, 637-653.
} 

\seealso{ \code{\link{MortSmooth_tpower}}. }

\examples{
x <- seq(0,1,length=500)
## B-splines basis of degree 1
B1 <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x),
                       ndx=10, deg=1)
matplot(x, B1, t="l", main="B-splines basis of degree 1")
## B-splines basis of degree 3
B3 <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x),
                       ndx=10, deg=3)
matplot(x, B3, t="l", main="B-splines basis of degree 3")
}

\keyword{smooth}
