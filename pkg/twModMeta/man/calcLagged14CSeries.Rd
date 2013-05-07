\name{calcLagged14CSeries}
\alias{calcLagged14CSeries}
\title{14C flux by applying lagged atmospheric 14C concentrations }
\description{Calculate the 14C flux by applying the atmospheric 14C concentration with a time lag.}
\usage{calcLagged14CSeries(inputYr, inputValue, delta14Catm, lag = 0, 
    iR14CStandard = c14Constants$iR14CStandard, delta13C = -25)}
\arguments{
  \item{inputYr}{the year of the carbon flux}
  \item{inputValue}{the value of the carbon flux}
  \item{delta14Catm}{dataframe columsn 1: yr and 2: delta14C signal of the atmosphere}
  \item{lag}{integer scalar: time lage in yr}
  \item{iR14CStandard}{isotopic ratio of the oxalic acid standard (heed that units of isotopes differ to scale to same magnitude)}
  \item{delta13C}{isotopic ratio of 13C in delta units for correction for fractionation}
}

\value{atomic ratio of input as multiple of the standard (dimensionless) corrected for fractionation.}

\author{Thomas Wutzler}





