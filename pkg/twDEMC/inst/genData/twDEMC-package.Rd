\name{twDEMC-package}
\alias{twDEMC-package}\alias{twDEMC}
\title{parallel DEMC}
\description{Parallel version of a DEMC differential evolution markov chain. DEMC is a variant of MCMC monte carlo markov chain.}
\details{The workhorse method for applying a Differential Evolution Markov Chain given a function of an unnormalized probability density is \code{\link{twDEMCBlockInt}}.It is usually invoked by \code{\link{twDEMCBlock.array}} or \code{\link{twDEMCBlock.twDEMCPops}}Further functionality of the twDEMC package deals with \itemize{\item{ Generating an initial population for twDEMC: \code{\link{initZtwDEMCNormal}}  }\item{ Transforming the results of twDEMC: \code{\link{subChains.twDEMC}}  }\item{ Transforming the parameter space: \code{\link{transOrigPopt.default}}  }\item{ Invoking fLogDen with proposal in a parallel load balanced way: \code{\link{twCalcLogDenPar}}  }\item{ Calculating marginal aggregates: \code{\link{marginals1d}}  }\item{ Laplace approximation of the normalizing constant: \code{\link{normConstLaplace}} and \code{\link{normConstLaplaceBridge}}  }\item{ Plotting: \code{\link{plotThinned.mcmc.list}}  }}%itemize}%details
\author{Thomas Wutzler}

\keyword{ package }


