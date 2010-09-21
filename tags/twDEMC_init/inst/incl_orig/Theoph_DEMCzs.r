rm(list=ls(all=TRUE))
library(nlme)
source("DEMC.ZS.r",echo= T)  # basic DE.MC sampler using joint updates, snooker and sampling differences from the past
source( "monitor.r", echo= FALSE)
# subject parameters together in groups of three (lKe, lKa, lCl)  rep(c(lKe, lKa, lCl), times =12)
# 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3
data(Theoph)

mySSfol <-function(x){ SSfol( x[1], x[2], x[3], x[4], x[5])}  # used in flogpost

Init.population.random <-function(N.subject,Npop = 1)  {
   X = NULL
   for (i in 1:Npop){
   # fixed parameters uniform between ranges
   lKe = runif(1,-2.95,-1.95)
   lKa = runif(1, -0.03,0.93)
   lCl = runif(1,-3.73,-2.73)
   # prior uniform on tau scale
   te =  runif(1,0.01 ,.1)
   ta =  runif(1,0.01, .1)
   tc =  runif(1,0.01, .1)
   ls2 =  runif(1,-1.19,-0.19)
   lKei = rnorm(N.subject, 0, te )
   lKai = rnorm(N.subject, 0, ta )
   lCli = rnorm(N.subject, 0, tc )
# subject parameters together in groups of three (lKe, lKa, lCl)  rep(lKe, lKa, lCl, times =12)
# 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3
   theta = c(lKe,lKa,lCl,(2*log(c(te,ta,tc))),ls2, as.vector(rbind(lKei, lKai,lCli)))
   X = cbind(X,theta)
   }
   rownames(X)= 1:43
   rownames(X)[1:7] =  c("lKe",  "lKa",  "lCl","lte2","lta2","ltc2", "ls2")
   X
}
                                   
#flogpost <-function(theta,MyData){
flogpost <-function(theta){        # cheaper to have MyData global
 # the function assumes each numbers of observations per subject (see X.tmp)
  N.subject=length(levels(MyData$Subject))   # 12
  NT = length(MyData$Dose)
  Npar = length(theta)
  N.time = NT / N.subject  # 11
  N.curve = 3   # number of parameters in curve  for an individual 
  Nfixed =  2*N.curve + 1
  #lKe = theta[1]
  #lKa = theta[2]
  #lCl = theta[3]
  #te = exp(theta[4]/2)  # tau _e
  #ta = exp(theta[5]/2)
  #tc = exp(theta[6]/2)
  sigma  = exp(theta[7]/2)
  sd.comp = exp(theta[(N.curve+1):(2*N.curve)]/2)        # sqrt of variance components
# subject parameters in Raneff together in groups of three (lKe, lKa, lCl)  rep(c(lKe, lKa, lCl), times =12)
# 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3    
  Raneff = theta[(Nfixed +1): Npar] 
  # parameters of indisvuals curves
  Ind.par = matrix(Raneff+  rep(theta[1:N.curve],N.subject), ncol= N.curve,byrow = T) # N.subject x N.curve
  #Ind.par.times=  matrix(apply(Ind.par, 1,rep, times = N.time),ncol=3,byrow = TRUE)
  X.tmp=  cbind(MyData$Dose,MyData$Time, matrix(apply(Ind.par, 1,rep, times = N.time),ncol=3,byrow = TRUE)) #Ind.par.times) 
  conc.expected  =  apply(X.tmp,1,mySSfol)
  loglikelihood  = sum(dnorm((MyData$conc - conc.expected), sd = sigma, log= TRUE))  
  # part of prior  from normally distributed random effectc
  log.pdf.raneff = sum(dnorm(Raneff,  sd= rep(sd.comp,N.subject), log= TRUE))    
 
  #  uniform prior for lKe, lKa, lCl
  # uniform prior on tau scale
  log.pdf.prior =  0.5* sum(theta[(N.curve+1):(2*N.curve)])
  loglikelihood + log.pdf.raneff + log.pdf.prior
 }
 
 

#ZmatL = as.matrix(read.table("Theoph_init.dat"))

Npar =  43
Npop = 3
MyData = Theoph
N.subject=length(levels(Theoph$Subject))   # 12
Z =   Init.population.random(N.subject, 10*Npar )
Nsim = (3*143333)/2    # halve of the number in the paper; takes ca. 80 minutes to run on 1.6 GHz processor
#Nsim = 3*14333
burnin = 0.2 * Nsim  
 
#Draws.ZS = DE_MC.ZS(Npop, Z, FUN=flogpost, X = Z[,1:Npop], pSnooker = 0.1, n.gen = floor(Nsim/Npop), n.burnin= floor(burnin/Npop), n.thin = max(floor(Nsim/100000),1)
#,eps.mult=0, eps.add = 0.001, MyData = Theoph)
Draws.ZS = DE_MC.ZS(Npop, Z, FUN=flogpost, X = Z[,1:Npop],pSnooker = 0.1, n.gen = floor(Nsim/Npop), n.burnin= floor(burnin/Npop), n.thin = max(floor(Nsim/100000),1)
,eps.mult=0, eps.add = 0)

dim(Draws.ZS$Draws) # 43 x ...
M = monitor.DE.MC(Draws.ZS$Draws[1:7,] , keep.all = TRUE , m= Npop)
  
M

proc.time()








