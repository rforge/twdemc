\name{modTwTwoDenEx1}
\alias{modTwTwoDenEx1}
\title{modTwTwoDenEx1}
\description{example model giving two predictions that can be compared to different observations}
\usage{modTwTwoDenEx1(theta, xSparse, xRich, thresholdCovar = 0)}
\arguments{
  \item{theta}{model parameters a and b}
  \item{xSparse}{numeric vector of Sparse input/output relationship}
  \item{xRich}{numeric vector of rich input/output relationship}
  \item{thresholdCovar}{model structural deficiency}
}
\details{model output y1 represents a longterm observations
It is based on longterm average of xRich instead of detailed values
\cr Model output y1 represents a short measurement campaing. 
During this campaing xSparse does not vary but detailed measurements of xRich are utilized 
\cr In the short-term relation, the model may simulate a detailed threshold in the covariate
or abstract from those details by thresholdCovar=0.}
\value{list with model predictions
\item{y1}{theta[1]*x1 + theta[2]*mean(xRich)/10
}}

\author{Thomas Wutzler}





