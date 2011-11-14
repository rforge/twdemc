DE.MC <- function(X, FUN, f = -2.38, n.iter = 10, n.thin = 1, n.burnin = 0, eps = 0, ...){
# Differential Evolution Markov Chain applied to X with logposterior specified by FUN
# X is the initial population: a matrix of number of parameters by number of individuals (k x Npop)
# f scaling parameter phi
# value
# FUN(theta ,...) with theta a k vector of parameter values and ... are other arguments to FUN (e.g. NULL, data, ..)
# Value
# $Draws k x (Npop*n) array  with n the number of retained simulations  [post process with monitor.DE.MC]
#         matrix(Draws[p,],nrow= Npop) is an Npop x n array (for each p in 1:k)
# $ accept.prob
# $ X.final
# $ logfitness.X.final
Npop = ncol(X)
Npar = nrow(X)
iseq = 1:Npop
F2 = abs(f)/sqrt(2*Npar)
if (f>0) F1 = F2 else F1 = 0.98		
accept = rep(NA,n.iter)
Draws = NULL
logfitness_X = apply (X, 2, FUN, ...)
for (iter in 1:n.iter) {
   accepti = 0
   if (iter%%10) F = F2 else F = F1		#longer jumps each ten occurences
   for (i in iseq){
     # select to random different individuals (and different from i) in rr, a 2-vector
     rr = sample(iseq[-i], 2, replace = FALSE)
     x_prop = X[,i] + F*(X[,rr[1]]-X[,rr[2]])  +  eps*rnorm(Npar,0,1)
     logfitness_x_prop = FUN(x_prop,  ...)
     if ((logfitness_x_prop - logfitness_X[i] )> log(runif(1)) ){
        accepti = accepti+1
        X[,i] = x_prop
        logfitness_X[i] = logfitness_x_prop
     }
  } # i loop
  accept[iter] = accepti
  if (!(iter%%n.thin) && iter > n.burnin) {  # retain sample
     Draws = cbind(Draws,X)
  }
} # n.iter
 list(Draws= Draws, accept.prob= accept/Npop, X.final = X, logfitness.X.final = logfitness_X)
}
