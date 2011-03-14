source(file.path(system.file(package="twKinrespTutorial"),"kinresp0.R"))

#for development of package only:
.tmp.f <- function(){
	library(nlme)
	source(file.path("R","kinrespModel1.R"))	# provided with package
	source(file.path("inst","kinresp0.R"))
}

#------------ least squares fit by hand (same as ML with assuming Gaussian errors)
#?optim		# help on the optim function
fCost1 <- function(
	### sum of sqared differences between y and pred
	par		##<< parameters: numeric vector (x0,r0,mumax)
	, ds	##<< dataset with columns resp and time
	, ...
){
	predResp <- kinrespModelVec( par, kinrespModel, time=ds$time, ...)
	# XX insert formula for calculating the log-Likelihood for Gaussian errors
	0
}
fitLS2 <- optim( par=c(x0=140, r0=0.003, mumax=0.24), fCost1, ds=ds )
fitLS2$par
coef(fitLS)		# compare to gnls fit

#take care: very senstive to local minima: 
# when starting from slightly different condtions, the otimization does not find the lower min
# need to check starting from different starting values (maybe grid)
fitLS3 <- optim( par=c(x0=140, r0=0.1, mumax=0.24), fCost1, ds=ds )
fCost1(fitLS3$par,ds)	# compare cost function with previous optimization
fCost1(fitLS2$par,ds)	 


#----------- MLE with assuming log-normal multiplicative error (constant relative error)
# see presentation for derivation: Yobs = Y*W		W~lognorm(0,sigma)
fCost2 <- function(
	### sum of sqared differences between y and pred
	par		##<< parameters: numeric vector (x0,r0,mumax)
	, ds	##<< dataset with columns resp and time
	, meanlog=-1/2*sdLog  
	, sdlog	  ##<< estimate of the standard deviation of the multiplicative error
){
	predResp <- kinrespModelVec( par, kinrespModel, time=ds$time)
	W <- ds$resp/predResp
	if( any(W<0) ) return( Inf )
	#XX: put formula here
}

# get an estimate of sd(W) by inspecting the gnls fit
W = ds$resp/fitted(fitLS)
#hist(W)		#maybe log-normally distributed
varW=var(W)
sdW=log(1+varW)
muW=-1/2*sdW
#library(debug)
#mtrace(cost2)
#fitLS4 <- optim( par=c(x0=140, r0=0.003, mumax=0.24), cost2, ds=ds, sd=sdLogW, method="BFGS" )
fitLS4 <- optim( par=c(x0=140, r0=0.003, mumax=0.24), fCost2, ds=ds, sdlog=sdW, meanlog=muW )
fitLS4$par			# compare the ML-lognorm-Hand estimates to the gnls estimates
coef(fitLS)
fitted4 <- kinrespModelVec( fitLS4$par, kinrespModel, time=ds$time)
plot( resp~time, data=ds)
lines(fitted(fitLS)~ds$time)
lines(fitted4~ds$time, col="red")

W2 <- ds$resp/fitted4
plot(W2 ~ log(fitted4))		# residuals 

#-------------- Profile Log-Likelihood
Lmin <- fCost1( coef(fitLS),ds=ds )
# get the search range from extending the confidence interval of gnls fit
tmp.cf <- tmp.cfExt <- confint(fitLS)	
tmp.cfr <- apply( tmp.cf, 1, diff )	
tmp.cfExt[,1] = tmp.cf[,1]-1.3*tmp.cfr 
tmp.cfExt[,2] = tmp.cf[,2]+1.3*tmp.cfr 

#----- ordinary profile holding r0 and x0 constant
# calculate the logLikelihood over a grid
# keep x0 and r0 constant, but vary mumax
nGrid=200
xrx <- expand.grid(x0 = coef(fitLS)["x0"]
	,r0 = coef(fitLS)["r0"]
	,mumax = seq(tmp.cfExt["mumax",1],tmp.cfExt["mumax",2],length.out=nGrid)
)
xrx$L <- apply(xrx[,c("x0","r0","mumax")], 1, fCost1, ds=ds ) 
plot( L ~ mumax, data=xrx)
# XX plot a line at the significance threshold

#two-dimensional case
xr <- expand.grid(x0 = seq(tmp.cfExt["x0",1],tmp.cfExt["x0",2],length.out=nGrid)
	,r0 = coef(fitLS)["r0"]
	,mumax = seq(tmp.cfExt["mumax",1],tmp.cfExt["mumax",2],length.out=nGrid)
)
xr$L <- apply(xr[,c("x0","r0","mumax")], 1, fCost1, ds=ds ) 

library(lattice)
#?levelplot
#levelplot( L ~ x0 + mumax, data=xr, col.regions=topo.colors(200))
levelplot( L ~ x0 + mumax, data=xr)

#only those values wihtin 95% cf, i.e. 
# XX select those cases within the 95% confidence range 
# xrc <- xr

#levelplot( L ~ x0 + mumax, data=xrc )
levelplot( L ~ x0 + mumax, data=xrc, col.regions=topo.colors(100) )
#levelplot( L ~ x0 + mumax, data=xrc, contour=TRUE)

#------ conditional profile: let r0 and x0 be optimized
fCost1Cond <- function(pars,ds){
	#pars <- xr[1,c("x0","r0","mumax")]
	# note the opitimization for x0 and r0 given mumax
	tmp <- optim( pars[c("x0","r0")], fCost1, ds=ds, parDefault=pars, method="BFGS" )
	tmp$value
}
#xr$Lc <- apply(xr[,c("x0","r0","mumax")], 1, fCost1Cond, ds=ds )
.tmp.f <- function(){
	# distribute calculation to 4 cores for speedup
	library(snowfall)
	sfInit(parallel=TRUE,cpus=4)
	sfExport("fCost1","kinrespModel","kinrespModelVec")
	xrx$Lc <- sfApply(xrx[,c("x0","r0","mumax")], 1, fCost1Cond, ds=ds )
}
xrx$Lc <- apply(xrx[,c("x0","r0","mumax")], 1, fCost1Cond, ds=ds )

plot( Lc ~ mumax, data=xrx, type="l", col="maroon")
lines( L ~ mumax, data=xrx, col="blue" )
abline(h=c(Lmin,Lmin+1.9), col="gray")

.tmp.f <- function(){
	#approximate log-Likelihood by kernel density estimate 	
	
	#library(snowfall)
	#sfInit(parallel=TRUE,cpus=4)
	#sfExport("fCost1","kinrespModel","kinrespModelVec")
	
	nGrid=50
	x <- expand.grid(x0 = seq(tmp.cfExt["x0",1],tmp.cfExt["x0",2],length.out=nGrid)
		,r0 = seq(tmp.cfExt["r0",1],tmp.cfExt["r0",2],length.out=nGrid)
		,mumax = seq(tmp.cfExt["mumax",1],tmp.cfExt["mumax",2],length.out=nGrid)
	)
	#x$L <- sfApply(x[,c("x0","r0","mumax")], 1, fCost1, ds=ds )
	x$L <- apply(x[,c("x0","r0","mumax")], 1, fCost1, ds=ds )
	
	#sfLibrary(np)
	#sfLibrary(twDEMC)
	library(np)
	library(twDEMC)
	form=L~x0+r0+mumax
	#mtrace(plot3DKDMarginals)
	tmp3 <- plot3DKDMarginals( form, x, dims=15, bandwidth=250, knotSpacing="equidistant", probs=c(0.25,0.5,0.75), nDrawPoints=500  )
	plot(tmp3, sample=NULL, nDrawPoints=250)
	plot(tmp3, nDrawPoints=800)
	plot(tmp3, probs=c(0.25,0.5,0.75,0.9), nDrawPoints=800)
	
	#tmp4 <- plot3DKDMarginals( form, sample0tvr, dims=10, bandwidth=250, nSample=1000, probs=c(0.5,0.75,0.95), nDrawPoints=250   )
	#plot( tmp4, probs=c(0.75,0.9,0.95), nDrawPoints=250   )
	#plot(tmp4,sample=NULL)
	#plot( tmp4, probs=c(0.5,0.75,0.95), nDrawPoints=250   )
	
	view3dTiltSpin(spin=60)	# and adjust zoom
	#rgl.bringtotop() 
	#rgl.snapshot(file.path("tmp","MarginalHowland_dO_h_tvrO.png"))
	play3dRound(12)
	#movie3dRound("MarginalHowland_dO_h_tvrO")
	
	open3d(windowRect=c(0,0,200,200)+20)	# adjust window widht
	plot(tmp3, probs=c(0.25,0.5,0.75,0.9), nDrawPoints=0, axes=FALSE)
	view3dTiltSpin(spin=60)	# and adjust zoom
	movie3dRound("Parameter_probability_Kinresp", 1/24, 6)
	
}



