\name{plot.Mort2Dsmooth}
\alias{plot.Mort2Dsmooth}

\title{ Plot of the Mort2Dsmooth outcomes }

\description{
  It produces two level plots, side-by-side, where the x- and y-axis are
  the predictors (\code{x} and \code{y}) of the given
  \code{Mort2Dsmooth} object and the numeric response is either the
  matrices of actual and fitted log-rates, or the matrices of actual and
  fitted counts. 
}

\usage{
\method{plot}{Mort2Dsmooth}(x,
                            type = c("logrates", "deaths"),
                            palette = c("rainbow",
                                        "heat.colors",
                                        "terrain.colors",
                                        "topo.colors",
                                        "cm.colors"), ...)
}

\arguments{
  \item{x}{ an object of class "Mort2Dsmooth", usually, a result of
  a call to \code{Mort2Dsmooth}. } 
  \item{type}{ the type of plot which should be returned. The
    alternatives are: \code{logrates} (default) and \code{death}. }
  \item{palette}{ a string with the color palette for creating a vector
  of contiguous colors. Default: \code{rainbow}. }
  \item{\dots}{ other plotting parameters passed to
  \code{levelplot}. Not in used. }
}

\details{
  The function provides a simple tool for checking the outcomes of a
  \code{Mort2Dsmooth} object. In a Poisson setting, one would commonly
  look at rates in log-scale, but actual and fitted counts can be
  plotted too. 
  
  The function uses \code{levelplot} from package \code{lattice} to
  construct the two level plots side-by-side. Plotting
  parameters are given and cannot be changed. The only expection is the
  color palette: five different palletes are available and described in
  the associated documentation: the default \code{rainbow} as
  well as \code{heat.colors}, \code{terrain.colors}, \code{topo.colors},
  and \code{cm.colors}.
  
  A palette with 9 colors is chosen and the breaks for the level
  plot are chosen according to the deciles of distributions of actual
  and fitted log-rates (counts).

}

\value{
  None. Function produces a plot.
}

\author{ Carlo G Camarda }

\seealso{ \code{\link{Mort2Dsmooth}} for computing
  \code{Mort2Dsmooth.object}. } 

\examples{
## selected data
ages <- 50:100
years <- 1970:2006
death <- selectHMDdata("Denmark", "Deaths", "Females",
                       ages = ages, years = years) 
exposure <- selectHMDdata("Denmark", "Exposures", "Females",
                          ages = ages, years = years)
## fit
fit <- Mort2Dsmooth(x=ages, y=years, Z=death,
                    offset=log(exposure),
                    method=3, lambdas=c(100,500)) 
## plotting log-death rates
plot(fit, palette="terrain.colors")
## plotting death counts
plot(fit, type="death")
}

\keyword{regression}
\keyword{smooth}
\keyword{hplot}
