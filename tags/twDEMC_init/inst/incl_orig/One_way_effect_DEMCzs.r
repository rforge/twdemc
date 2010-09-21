#
# ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version of
#    the genetic algorithm Differential Evolution: easy Bayesian computing
#    for real parameter spaces. Statistics and Computing, 16, 239-249.
#
# and
#
# ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain 
#  with snooker updater and fewer chains. Statistics and Computing
#    http://dx.doi.org/10.1007/s11222-008-9104-9 .
#
# Simple two-dimensional example of DE-MC based on section 4.1 (one-way random-effects model) in both papers
# Here we use the marginal pdf [function flogpost] of log(sigma^2) and log(tau^2) to sample from
# In the paper the problem was 7-dimensional as the pdf  included mu and the 4 theta's
# Hence d = 2 
# Note: for low dimensional examples there seems little advantage of DE.MC.ZS compared to DE.MC

# Functions that the user must specify are
#  the logposterior , here the function flogpost
#  the initial population, here the function Init.population

rm(list=ls(all=TRUE))

source("DEMC.r")    # basic DE.MC sampler using joint updates
source("DEMC.ZS.r",echo= T)  # basic DE.MC sampler using joint updates, snooker and sampling differences from the past

source("monitor.r")  #  postprocessor monitor.DE.MC adapted from Adrew Gelman's monitor

original.data=matrix(c(
 0.34,0.12,1.23,0.70,1.75,0.12,
0.91,2.94,2.14,2.36,2.86,4.55,
6.31,8.37,9.75,6.09,9.82,7.24 ,
17.15,11.82,10.95,17.20,14.35,16.82),
  nrow= 4, ncol= 6, byrow =TRUE)

# Liu & Hodges (2003) say that they analyze the sqrt of original data matrix
ymat = sqrt(original.data)

theta = c(10,20,30,1,2,3,4)

# Functions that the user must specify
flogpost.d7 <- function(theta, MyData, alpha= 1,beta=10,a=1.85,b=0.1, I = 6, J = 4){
# logposterior of 7 parameters mu, log(sigma^2)  and log(tau^2) , theta1,...,theta4
# data(J,I)
I = 6  # number of ind per group
#J = 4  # number of groups
# for single values of logs2
  mu = theta[1]
  logs2 = theta[2]
  logt2 = theta[3]
  s2 = exp(logs2)
  t2 = exp(logt2)
  raneff = theta[4:7]
  loglikelihood = sum(dnorm( MyData,mean = matrix(rep(raneff+ mu,I),ncol= I),sd = sqrt(s2), log = TRUE))  
  logprior.raneff =   sum(dnorm(raneff,sd = sqrt(t2), log = TRUE))  
  minlogprior.vcomp =  alpha * logs2 + beta/s2 + a * logt2 +  b/t2
  loglikelihood + logprior.raneff - minlogprior.vcomp
}


Init.population.d7 <- function(Npop, alpha= 1,beta=10,a=1.85,b=0.1, I = 6, J = 4){ 
# the population is a matrix of #par  by #individual, here 2x Npop
# Here the initial population is drawn from the prior
 mu = runif(Npop, min = -5,max = 5)
 sigma2.inv = rgamma(Npop, alpha, rate = beta)
 tau2.inv = rgamma(Npop, a, rate = b)
 tau = sqrt(1/tau2.inv)
 theta1_4 = matrix(rnorm(J*Npop,mu,sd = tau), ncol= Npop)
 rbind(mu,-log(sigma2.inv), -log(tau2.inv),theta1_4)
}

Init.population.d2 <- function(Npop, alpha= 1,beta=10,a=1.85,b=0.1, I = 6, J = 4){ 
# the population is a matrix of #par  by #individual, here 2x Npop
# Here the initial population is drawn from the prior
 sigma2.inv = rgamma(Npop, alpha, rate = beta)
 tau2.inv = rgamma(Npop, a, rate = b)
 rbind(-log(sigma2.inv), -log(tau2.inv))
}


# Functions that the user must specify
flogpost.d2 <- function(theta, data, alpha= 1,beta=10,a=1.85,b=0.1, I = 6, J = 4){
# logposterior of log(sigma^2)  and log(tau^2) [after integrating analytically over mu and theta]
# data(J,I)
#I = 6  # number of ind per group
#J = 4  # number of groups
# for single values of logs2
  logs2 = theta[1]       
  logt2 = theta[2]
  s2 = exp(logs2)
  t2 = exp(logt2)
  diag_elem = ((I-1)*t2+s2)/(s2*(I*t2 + s2))
  offdiag = -t2/(s2*(I*t2 + s2))
  invCII =  matrix(offdiag,nrow =I,ncol =I)
  diag(invCII) =  diag_elem

  ydif = data - mean(data)
  detCII = I*t2*s2^(I-1) + s2^I  # det(C)
  Vmu = (s2/I + t2)/J
  logpdf =   0.5* (log(Vmu) -J*log(detCII) - sum(diag(ydif%*%invCII%*% t(ydif)))   )

  minlogprior =  alpha * logs2 + beta/s2 + a * logt2 +  b/t2
  logpdfst = logpdf - minlogprior
  logpdfst
}






Nsim =  10000
burnin = 0.2 *Nsim
# two level normal random effects model

# original DE.MC method (ter Braak, 2006, Stat.Comp. 16:239-249, http://dx.doi.org/10.1007/s11222-006-8769-1)
d = 2  # number of parameters (length of theta) in flogpost
Npop = 3*d
X = Init.population.d7(Npop)

Draws = DE.MC(X, flogpost.d7, n.iter = floor(Nsim/Npop),
                     n.burnin= floor(burnin/Npop),MyData = ymat)
monitor.DE.MC(Draws$Draws, keep.all = TRUE , m= Npop)

logksi = Draws$Draws[2,] -   Draws$Draws[3,]
phi = exp(Draws$Draws[2,]) /(6*exp(Draws$Draws[3,])+  exp(Draws$Draws[2,]))
Plogksi = quantile(logksi, prob= c(.025, .5, .975))
Pphi = quantile(phi, prob= c(.025, .5, .975))


# DE.MC.ZS method (ter Braak & Vrugt, 2008, Stat.Comp. http://dx.doi.org/10.1007/s11222-008-9104-9 )
d = 7   # number of parameters (length of theta) in flogpost
Npop = 3
m0 = 5*d     # the default of 10d is worse for low d
Z =  Init.population.d7(m0)
Draws.ZS = DE_MC.ZS(Npop, Z, FUN=flogpost.d7, X=X[,1:Npop]
	, pSnooker = 0.1, n.gen = floor(Nsim/Npop), n.burnin= floor(burnin/Npop), n.thin = max(floor(Nsim/100000),1)
	,eps.mult=0, eps.add = 0.001, MyData = ymat)

monitor.DE.MC(Draws.ZS$Draws, keep.all = TRUE , m= Npop)

logksi.new = Draws.ZS$Draws[2,] -   Draws.ZS$Draws[3,]
phi.new = exp(Draws.ZS$Draws[2,]) /(6*exp(Draws.ZS$Draws[3,])+  exp(Draws.ZS$Draws[2,]))

Plogksi.new = quantile(logksi.new, prob= c(.025, .5, .975))
Pphi.new = quantile(phi.new, prob= c(.025, .5, .975))

results = rbind(Plogksi,Plogksi.new,Pphi,Pphi.new)
results


