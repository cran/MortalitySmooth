\name{MortSmooth_tpower}

\alias{MortSmooth_tpower}

\title{Truncated p-th Power Function}

\description{
  This is an internal function of package \code{MortalitySmooth} which
  constructs a truncated \emph{p}-th power function along an abscissa
  within the function \code{MortSmooth_bbase} 
}

\usage{
MortSmooth_tpower(x, t, p)
}

\arguments{
  \item{x}{ vector for the abscissa of data. }
  \item{t}{ vector of truncation points. }
  \item{p}{ degree of the power. }
}

\details{
  Internal function used in \cite{\code{MortSmooth_bbase}}. The vector
  \code{t} contains the knots. The simplest system of truncated power
  functions uses \code{p = 0} and it consists of step functions with
  jumps of size 1 at the truncation points \code{t}. 
}

\references{ Eilers and Marx (2004). Splines, Knots, and
  Penalties. Unpublished manuscript. } 

\author{ Carlo G Camarda }

\seealso{ \code{\link{MortSmooth_bbase}}. }

\examples{
x <- seq(0,1,length=100)
f1 <- MortSmooth_tpower(x=x, t=0.1, p=3)
f2 <- MortSmooth_tpower(x=x, t=0.2, p=3)
f3 <- MortSmooth_tpower(x=x, t=0.3, p=3)
## Truncated 3rd power functions equally-spaced
plot(x, f1, t="l",
     main="Truncated 3rd power functions equally-spaced") 
lines(x, f2, col=2)
lines(x, f3, col=3)
}

\keyword{smooth}
