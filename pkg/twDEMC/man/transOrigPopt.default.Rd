\name{transOrigPopt.default}
\alias{transOrigPopt.default}
\title{transOrigPopt default}
\description{Transform vectors from normal to original scale.}
\usage{
\method{transOrigPopt}{default}(normpopt, poptDistr = eval(parse(text = "parDistr$trans[names(normpopt)]")), ...)
}
\arguments{
  \item{normpopt}{numerical vector/array with values at transformed, i.e. normal, scale}
  \item{poptDistr}{character vector/array of kind of transformation ("lognorm"/"logitnorm")
values with other characters indicate no transformation
default assumes vector parDistr$trans in environement }
  \item{\dots}{
}
}
\details{This generic method is provided for in several forms for first argument. \itemize{
\item{ as a named vector: this method  } 
\item{ as a matrix: \code{\link{transOrigPopt.matrix}}  } 
\item{ as a coda's mcmc.list: \code{\link{transOrigPopt.mcmc.list}}  }
\item{ as a list of class twDEMC \code{\link{transOrigPopt.twDEMC}}  }
}

There are further methods deal with parameter transformations. \itemize{
\item{ transforming from original to normal scale: \code{\link{transNormPopt.default}}  } 
\item{ calculating mu and sigma at normal scale from quantiles: \code{\link{twQuantiles2Coef}}  } 
\item{ constraing the result list of \code{\link{twQuantiles2Coef}}, and adding variance-covariance matrix: \code{\link{twConstrainPoptDistr}}  }
}

Argument \code{poptDistr} should have the same dimensions as normpopt. However, it is recycled
By this way it is possible to specify only one value, or vector corresponding to the rows of a matrix.}
\value{Normpopt with some values transformed by exp (poptDist=="lognorm") or plogis (poptDistr=="logitnorm").}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}
\code{\link{transNormPopt.default}}}

