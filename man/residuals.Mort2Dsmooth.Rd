\name{residuals.Mort2Dsmooth}
\alias{residuals.Mort2Dsmooth}

\title{ Extract 2D P-splines Model Residuals }

\description{
  Extracting different types of residuals from a \code{Mort2Dsmooth}
  object. 
}

\usage{
\method{residuals}{Mort2Dsmooth}(object,
                                 type = c("deviance", "pearson", "anscombe", "working"), ...)
}

\arguments{
  \item{object}{ An object of class "Mort2Dsmooth", usually, a result of
  a call to \code{Mort2Dsmooth}. } 
  \item{type}{ The type of residuals which should be returned. The
  alternatives are: "deviance" (default), "anscombe", "pearson" and
  "working". } 
  \item{\dots}{ Further arguments passed to or from other methods. } 
}

\details{
  The references define the types of residuals. 
  
  The way of computing the residuals are described in Section 2.4 of
  McCullagh and Nelder's book. The working residuals are merely the
  differences between fitted and actual counts. 
}

\value{
  A matrix of the selected type of residuals over both the \code{x} and
  the \code{y} axes in the \code{Mort2Dsmooth} object. 
}

\references{ 
  Davison, A. C. and Snell, E. J. (1991). \emph{Residuals and
  diagnostics.}  In: Statistical Theory and Modelling. In Honour of Sir
  David Cox, FRS, eds. Hinkley, D. V., Reid, N. and Snell, E. J.,
  Chapman \& Hall. 

  McCullagh P. and Nelder, J. A. (1989). \emph{Generalized Linear
 Models.} London: Chapman and Hall.
 }

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}} for computing
  \code{Mort2Dsmooth.object}, \code{\link{residuals}}. } 

\examples{
# selected data
ages <- 30:80
years <- 1970:2006
death <- selectHMDdata("Switzerland", "Deaths", "Males",
                       ages = ages, years = years) 
exposure <- selectHMDdata("Switzerland", "Exposures", "Males",
                          ages = ages, years = years)
# fit
fit <- Mort2Dsmooth(x=ages, y=years, Z=death, offset=log(exposure),
                    method=3, lambdas=c(300,10))

# extracting residuals
devR <- resid(fit, type="deviance")
ansR <- resid(fit, type="anscombe")
peaR <- resid(fit, type="pearson")
worR <- resid(fit, type="working")

# plotting deviance residuals over age and years
res.list <- list(ages=ages, years=years)
res.grid <- expand.grid(res.list)
res.grid$dev <- c(devR)
levelplot(dev~years*ages, res.grid,
          at=c(min(devR), -2, -1, 1, 2, max(devR)))
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
