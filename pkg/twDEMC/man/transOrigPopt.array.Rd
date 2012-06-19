\name{transOrigPopt.array}
\alias{transOrigPopt.array}
\title{transOrigPopt array}
\description{Applies \code{\link{transOrigPopt.default}} to each column, i.e. variable, of \code{normopt}.}
\usage{
\method{transOrigPopt}{array}(normpopt,
poptDistr = eval(parse(text = "parDistr$trans[rownames(normpopt)]")),
...)
}
\arguments{
  \item{normpopt}{numerical matrx with values at transformed, i.e. normal, scale}
  \item{poptDistr}{character vector of kind of transformation ("lognorm"/"logitnorm") for each column of normpopt}
  \item{\dots}{
}
}
\details{either poptDistr has names for each column name
or poptDistr has the same length as colnames(normpopt)}


\author{Thomas Wutzler}





