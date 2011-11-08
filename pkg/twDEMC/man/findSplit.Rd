\name{findSplit}
\alias{findSplit}
\title{findSplit}
\description{determine the parameter and its value, where to split the parameter space}
\usage{
findSplit(aSample, nSplit = 4, rVarCrit = 10^2, rAlphaSlopeCrit = base:::pi/3, pSlopeMax = 0.05, iVars = 1:ncol(aSample), jVars = 1:ncol(aSample), isBreakEarly = TRUE)
}
\arguments{
  \item{aSample}{the sample (stacked chains)}
  \item{nSplit}{how many points per parameter dimension to check for split}
  \item{rVarCrit}{minimla ratio of variances of a variable between right and left side of a splitting point}
  \item{rAlphaSlopeCrit}{minimal ratio of angle between scaled slopes, defaults to a third of a half cirle}
  \item{pSlopeMax}{confidence level for the slope, above which no difference in slopes is calculated}
  \item{iVars}{integer vector: index or parameter dimensions, i.e. variables, to check for split}
  \item{jVars}{integer vector: index or parameter dimensions, i.e. variables, to check for different scales of variance}
  \item{isBreakEarly}{if TRUE returns on first checked variable that has a breakpoint, else check all variables and return order of variables}
}
\details{First it checks for different scales of variance in other variables.

Next it checks for different angles of the normalized slopes in relation of the parameters }
\value{list with components  
\item{split}{named scalar: splitting value with the name of the parameter dimension that is to split}
\item{varName}{name of the splitting variable}
\item{perc}{scalar: percentile of the splitting point}
\item{iVars}{integer vector: argument \code{iVars} ordered by maximum proportion of variances or difference in slopes}
\item{jVars}{integer vector: argument \code{jVars} ordered by proportion of varainces or difference in slopes for iVars[1]}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(den2dCorTwDEMC)
ss1 <- ss <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
#mtrace(findSplit)
(res <- res0 <- findSplit(ss1, isBreakEarly=FALSE))	# check all variable and return ordering
(res <- res1 <- findSplit(ss1))

# successively find splits in subsets
ss2 <- ss <- ss1[ ss1[,res$varName] > res$split,]
res2 <- findSplit(ss2)
ss3 <- ss <- ss2[ ss2[,res2$varName] > res2$split,]
res3 <- findSplit(ss3)
# NA indicates: no further split found

# try different orders of the variables
#mtrace(findSplit)
(res <- findSplit(ss1, iVars=c(2,1) ))
(res <- findSplit(ss1, iVars=c("b") ))
plot( ss1[,1], ss1[,2], col=c("blue","red")[as.integer(ss1[,res1$varName] >= res1$split)+1])
}
