\name{selectHMDdata}
\alias{selectHMDdata}

\title{Selecting Data from the HMDdata Object}

\description{
Creates subset of the HMDdata object.
}

\usage{
selectHMDdata(country,
              data = c("Population", "Deaths", "Exposures", "Rates"), 
              sex = c("Females", "Males", "Total"), 
              ages = NULL, years = NULL)
}

\arguments{
  \item{country}{ character string for the country name. }
  \item{data}{ character string showing type of demographic series:
  either "Population", "Deaths", "Exposures" or "Rates". Default value
  "Population". } 
  \item{sex}{ character string showing sex of demographic series: either
  "Females", "Males" or "Total". Default value "Females". } 
  \item{ages}{ vector of ages to extract from data. If it is
  \code{NULL}, all ages will be selected. } 
  \item{years}{ vector of years to extract from data. If it is
  \code{NULL}, all available years will be selected. } 
}

\details{
  \code{HMDdata} are loaded with the package.
  
  Country names is not case-sensitive and initials can be used. However,
  distinguishable country name is necessary. Type of \code{data} and
  \code{sex} are case-sensitive and initials can be used. 
  
  Available ages are always from 0 to 110, whereas years depends on the
  country. Warning massages will be provided in case of selected default
  values. Attributes about country, data and sex will be given to the
  outcome. 
}

\value{
  Matrix object from \code{\link{HMDdata}} with a subset of country,
  type of data, sex, ages and years. 
}

\author{Carlo G Camarda}

\seealso{ \code{\link{HMDdata}}. }

\examples{
## "Complete" example
x <- 1950:2000
den50 <- selectHMDdata("Denmark", "Death", "Females", 50, x)
plot(x, den50, main="Danish female deaths at age 50")
## "Incomplete" example with warning
jap50 <- selectHMDdata("jap", "Pop", "F", 50)
## Wrong example (don't run)
## selectHMDdata("Sw", "Pop", "F", 50, 2000)
## Sw can stand for both Sweden and Switzerland
}

\keyword{datasets}
