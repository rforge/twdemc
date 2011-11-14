# Example R code for DE_MC.ZS
#
#  ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain 
#  with snooker updater and fewer chains. Statistics and Computing
#    http://dx.doi.org/10.1007/s11222-008-9104-9 .
#
# The multiplicative error is detailed in: 
#  Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       
#       Accelerating Markov chain Monte Carlo simulation by differential evolution with         
#       self-adaptive randomized subspace sampling, International Journal of Nonlinear          
#       Sciences and Numerical Simulation, In Press. 2008.  
#
# see also  
# ter Braak, C. J. F. (2006). A Markov Chain Monte Carlo version of
#    the genetic algorithm Differential Evolution: easy Bayesian computing
#    for real parameter spaces. Statistics and Computing, 16, 239-249.
#


files 

Function files:
DEMC.r                    - simplest implementation of the basic DE-MC sampler from Ter Braak 2006 (without snooker, without sampling difference from the past
DEMC.ZS.r                 - implementation of ter Braak & Vrugt 2008 DE-MC with snooker and sampling differences from the past
monitor.r                 - functions adapted from Andrew Gelman's R2Bugs for monitoring the output of the sampler



Example files 
One_way_effect_DEMCzs.r   - example section 4.1     of ter Braak & Vrugt 2008
Theoph_DEMCzs.r           - example section 4.2     of ter Braak & Vrugt 2008



