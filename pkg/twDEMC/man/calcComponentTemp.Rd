\name{calcComponentTemp}
\alias{calcComponentTemp}
\title{calcComponentTemp}
\description{calculating the temperature of result components of logDensity}
\usage{calcComponentTemp(temp,
 TFix, TProp, useMultiT = TRUE, posTFix = match(names(TFix), 
    names(TProp)))}\arguments{
  \item{temp}{numeric scalar >= 1: global temperature}
  \item{TFix}{named vector: temperature for the components that does not change but is held fixed}
  \item{TProp}{named numeric vector [0,1]: temperature proportions of result components determines names and lenght of result}
  \item{useMultiT}{if set to FALSE only only global temperatue and TFix are applied}
  \item{posTFix}{position index of TFix in vector of temperatures, specifiy for performance}
}



\author{Thomas Wutzler}





