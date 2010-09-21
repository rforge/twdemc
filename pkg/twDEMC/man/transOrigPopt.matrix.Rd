\name{transOrigPopt.matrix}
\alias{transOrigPopt.matrix}
\title{transOrigPopt matrix}
\description{Applies \code{\link{transOrigPopt.default}} to each column of \code{normopt}.}
\usage{
transOrigPopt.matrix(normpopt, poptDistr = eval(parse(text = "parDistr$trans[colnames(normpopt)]")), ...)
}
\arguments{
  \item{normpopt}{numerical matrx with values at transformed, i.e. normal, scale}
  \item{poptDistr}{character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt}
  \item{\dots}{
}
}
\details{either poptDistr has names for each column name
or poptDistr has the same length as colnames(normpopt)}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





