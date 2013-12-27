\name{residuals.Mort1Dsmooth}
\alias{residuals.Mort1Dsmooth}

\title{ Extract 1D P-splines Model Residuals }

\description{
  Extracting different types of residuals from a \code{Mort1Dsmooth}
  object. 
}

\usage{
\method{residuals}{Mort1Dsmooth}(object,
          type = c("deviance", "pearson",
                   "anscombe", "working"), ...)
}

\arguments{
  \item{object}{ an object of class "Mort1Dsmooth", usually, a result of
  a call to \code{Mort1Dsmooth}. } 
  \item{type}{ the type of residuals which should be returned. The
  alternatives are: "deviance" (default), "anscombe", "pearson" and
  "working". }
  \item{\dots}{ Further arguments passed to or from other methods. Not
  in used. Not in use. } 
}

\details{
  The references define the types of residuals.
  
  The way of computing the residuals are described in Section 2.4 of
  McCullagh and Nelder's book. The working residuals are merely the
  differences between fitted and actual counts. 
}

\value{
  A vector of the selected type of residuals for each of the predictor
  in the \code{Mort1Dsmooth} object. 
}

\author{ Carlo G Camarda }

\references{ 
  Davison, A. C. and Snell, E. J. (1991). \emph{Residuals and
  diagnostics.}  In: Statistical Theory and Modelling. In Honour of Sir
  David Cox, FRS, eds. Hinkley, D. V., Reid, N. and Snell, E. J.,
  Chapman & Hall. 

  McCullagh P. and Nelder, J. A. (1989). \emph{Generalized Linear
  Models.} London: Chapman & Hall.
}

\seealso{ \code{\link{Mort1Dsmooth}} for computing
  \code{Mort1Dsmooth.object}. } 

\examples{
## selected data
years <- 1970:2006
death <- selectHMDdata("Denmark", "Deaths", "Females",
                       ages = 60, years = years)
exposure <- selectHMDdata("Denmark", "Exposures", "Females",
                          ages = 60, years = years)
## fit
fit <- Mort1Dsmooth(x=years, y=death, offset=log(exposure),
                    method=3, lambda=1000)

## extracting residuals
devR <- resid(fit, type="deviance")
ansR <- resid(fit, type="anscombe")
peaR <- resid(fit, type="pearson")
worR <- resid(fit, type="working")
## summaries
summary(devR)
summary(ansR)
summary(peaR)
summary(worR)
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
