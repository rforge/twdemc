\name{squeeze.twDEMCPops}
\alias{squeeze.twDEMCPops}
\title{squeeze twDEMCPops}
\description{Reduces the rows of populations so that all populations have the same number of samples. }
\usage{\method{squeeze}{twDEMCPops}(x, ..., length.out = min(getNSamples(x)), 
    dropShortPops = FALSE)}
\arguments{
  \item{x}{the twDEMCPops list to thin}
  \item{\dots}{
}
  \item{length.out}{number vector (nPops) of steps in each population, or numeric scalar specifying the lenght for all populations}
  \item{dropShortPops}{if set to TRUE, pops with less samples than length.out[1] are dropped}
}
\details{all populations with 

components \code{thin,Y,nGenBurnin} are kept, but may be meaningless after subsetting.}


\author{Thomas Wutzler}



\seealso{\code{\link{subset.twDEMCPops}}}
\examples{
if( FALSE){
	tmp <- squeeze(res)
	getNSamples(tmp)
	tmp2 <- subPops(res,2)
	#mtrace(squeeze.twDEMCPops)
	getNSamples( squeeze(tmp2,length.out=10) )
	getNSamples( squeeze(tmp,length.out=10) )	# should produce error
	getNSamples( squeeze(tmp,length.out=10, dropShortPops=TRUE) )	# should produce error
} 
}
