\name{findSplit}
\alias{findSplit}
\title{findSplit}
\description{determine the parameter and its value, where to split the parameter space}
\usage{
findSplit(aSample, nSplit = 4, rVarCrit = 3^2, rAlphaSlopeCrit = base:::pi/4, pSlopeMax = 0.05, iVars = 1:ncol(aSample), jVarsVar = 1:ncol(aSample), jVarsSlope = 1:ncol(aSample), isBreakEarly = TRUE, checkSlopesFirst = data.frame(ivar = 1, j1AlphaSlope = 1, jAlphaSlope = 2)[FALSE, ], maxNSample = 128, calcVarParallel = TRUE, debugSequential = FALSE)
}
\arguments{
  \item{aSample}{the sample (stacked chains)}
  \item{nSplit}{how many points per parameter dimension to check for split}
  \item{rVarCrit}{minimal ratio of variances of a variable between right and left side of a splitting point}
  \item{rAlphaSlopeCrit}{minimal ratio of angle between scaled slopes, defaults to a third of a half cirle}
  \item{pSlopeMax}{confidence level for the slope, above which no difference in slopes is calculated}
  \item{iVars}{integer vector: index or parameter dimensions, i.e. variables, to check for split}
  \item{jVarsVar}{integer vector: index or parameter dimensions, i.e. variables, to check for different scales of variance}
  \item{jVarsSlope}{integer vector: index or parameter dimensions, i.e. variables, to check for differences in angle of normalized slopes}
  \item{isBreakEarly}{if TRUE and argument \code{checkSlopesFirst} is given, then check slope angles for variables given in checkSlopesFirst first and if break is found, do not evaluate all the other slope angles}
  \item{checkSlopesFirst}{data.frame with entries ivar, j1AlphaSlope and j2AlphaSlope as returned in entry resD in result of ressplit}
  \item{maxNSample}{if given a value, then aSample is thinned before to given number of records (for efficiently calculating variances)}
  \item{calcVarParallel}{having less than 30 parameters it is often faster not to parallelize variance calculation}
  \item{debugSequential}{by default calcualation is distributed to sup-processes for parralel execution, set to TRUE for non-distributed execution}
}
\details{First it checks for different scales of variance in other variables.

Next it checks for different angles of the normalized slopes in relation of the parameters}
\value{list with components  
\item{split}{named scalar: splitting value with the name of the parameter dimension that is to split}
\item{varName}{name of the splitting variable}
\item{perc}{scalar: percentile of the splitting point}
\item{resD}{dataframe of result details giving for each split variable the indices of variable with maximum proportion of variances and and the variable with maximum angle in correlation}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
#twUtestF(findSplit) # there are unit tests for this function
data(den2dCorTwDEMC)
ss1 <- ss <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
#mtrace(findSplit)
(res <- res0 <- findSplit(ss1))	# returns before checking slope angles
(res <- res1 <- findSplit(ss1, rVarCrit=Inf))	# find slope angles

# successively find splits in subsets
ss2 <- ss <- ss1[ ss1[,res$varName] > res$split,]
res2 <- findSplit(ss2)
ss3 <- ss <- ss2[ ss2[,res2$varName] > res2$split,]
res3 <- findSplit(ss3)

# try different orders of the variables
(res <- findSplit(ss1, iVars=c("b") )) # NA indicates: no further split found

#visualize the split
plot( ss1[,"a"], ss1[,"b"], col=c("blue","red")[as.integer(ss1[,res0$varName] >= res0$split)+1])
}
