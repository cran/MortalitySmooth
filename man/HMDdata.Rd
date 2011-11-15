\name{HMDdata}

\alias{HMDdata}

\docType{data}

\title{Population and Mortality Data}

\description{
  Age-specific population, deaths, exposures and rates from the Human
  Mortality Database. 
}

\usage{data(HMDdata)}

\format{Object of class \code{HMDdata} containing, for each country, the
  following components: 

\describe{
\item{country}{ country name. }
\item{year}{ vector of years. }
\item{age}{ vector of ages. }
\item{pop}{ list of matrices containing population with one age group
  per row and one column per year. Matrices: \code{female}, \code{male},
  \code{total}. } 
\item{death}{ list of matrices containing deaths in same form as
  \code{pop}. } 
\item{exposure}{ list of matrices containing exposure in same form as
  \code{pop}. } 
\item{rate}{ list of matrices containing mortality rate in same form as
  \code{pop}. } 
}}

\details{
  The list contained 4 countries (Denmark, Japan, Sweden and
  Switzerland). Data taken from the Human Mortality Database on 30 August
  2011. Data are accessible either manually or using
  \code{\link{selectHMDdata}}.
  
  For a given country, the matrices have the same dimensions. The vector
  \code{age} is the same for each country (0:110), whereas years depends
  on the availability. 
  
  \code{pop} indicates population size on January 1st of each age and
  year. 
  
  \code{death} are death counts occurred during one age-year interval. 
  
  \code{exposure} population is approximated by the average of the
  population size in the beginning and at the end of the year. 

  \code{rate} are computed as ratio between death counts and exposures. 
}

\author{Carlo G Camarda}

\seealso{ \code{\link{selectHMDdata}}. }

\source{Human Mortality Database \url{www.mortality.org}}

\references{Human Mortality Database (2011). University of California,
  Berkeley (USA), and Max Planck Institute for Demographic Research
  (Germany). Available at \url{www.mortality.org}.} 

\examples{
data(HMDdata)
plot(HMDdata$Japan$year, HMDdata$Japan$pop$female[81,],
     main="Japan population, female, age 80") 
## list of available countries
names(HMDdata)
}

\keyword{datasets}
