tmp.f <- function(){
	library(abind)
	library(coda)
	library(snowfall)
	sfInit(parallel=TRUE, cpus=4)
	#sfInit(parallel=FALSE)
	#sfStop()
}




#psrs <- list( psrb_s2FGS, psrb_s2FGS )



#InitDEMCzsp <- function(thetaPrior, covarTheta, nPop=4, m0=max(4,(8*length(tthetaPrior)%/%Npop) ){






tmp.f <- function(){
	#testing temperature adaptation calculation
	pAcceptWindowWidth = 100
	T0=100
	T0=1000
	T0=10
	ctrl <- list(thin=8)
	TcurStep=T0
	
	tmp.n <- 8*pAcceptWindowWidth
	tmp.T0 <- max(10,T0)
	Tb <- tmp.T0
	Ta <- (log(Tb)-1)/tmp.n

	#Ta=log(/1)/(4*pAcceptWindowWidth-1); Texpa <- exp(Ta)	#max needed to give it a chance to increase Temp when starting from 1
	#Tb=max(10,T0)*Texpa
	
	iT <- (log(Tb/TcurStep)/Ta)		# position on temperature curve
	Tb*exp( -0.66*Ta*(iT-(1:ctrl$thin)) )	
	tmp <- (-tmp.n:tmp.n); plot(Tb*exp( -Ta*tmp)~tmp); abline(v=0); abline(h=TcurStep)
	tmp <- (-10:10); plot(Tb*exp( -Ta*tmp)~tmp); abline(v=0); abline(h=TcurStep)
}


tmp.f <- function(){
	#former part of DEMCzsp
	# if acceptance rate is too low, postpone cooling
	# !caution: at the end Tend will not be reached
	tmp <- iGen + (1:ctrl$thin)
	if( median(pAccept[mZ,]) < 0.1 ){
		Tstep[(iGen+1):nGen] <- c( rep(Tstep[iGen],ctrl$thin), (if(iGen<nGen) Tstep[(iGen+1):(nGen-ctrl$thin)] else c()) )
	}
	
	# former part of DEMCzsp
	# only complicates it, do not use call argument any more
	if( !is.null(Zinit$call) ){	#update parameters from former call
		fcl <- Zinit$call
		if( !hasArg(fLogLik) ){
			ccl$fLogLik <- fcl$fLogLik
			fLogLik <- eval.parent(fcl$fLogLik)
		}
		if( !is.null(fcl$DEMCzspControl) ){
			#reuse old control arguments that have not been respecified
			fCtrl <- eval.parent(fcl$DEMCzspControl)
			argNames <- names(fCtrl)[ !(names(fCtrl) %in% names(DEMCzspControl)) ]
			ctrl[argNames] <- fCtrl[argNames]
			# do not make changes permannet
			#ccl$DEMCzspControl <- as.call( c(list, as.list(ccl$DEMCzspControl)[-1], as.list(fcl$DEMCzspControl)[argNames]) )
			ccl$DEMCzspControl <- fcl$DEMCzspControl
		}
		if( !is.null(fcl$argsFLogLik) ){
			#reuse old arguments to fLogLik that have not been respecified, when specified as a name, evaluate first
			fArgs <- eval.parent(fcl$argsFLogLik)
			if( is.list(fArgs) ){
				argNames <- names(fArgs)[ !(names(fArgs) %in% names(argsFLogLik)) ]
				argsFLogLik[argNames] <- fArgs[argNames]
			}else{
				if( !hasArg(argsFLogLik) )
					argsFLogLik <- fArgs	
			}
			#ccl$argsFLogLik <- substitute(argsFLogLik)	#will replace reference to list by contents of list
			# do not update, provide the former (and by recursion first argument set)
			#ccl$argsFLogLik <- as.call( c(list, as.list(ccl$argsFLogLik)[-1], as.list(fcl$argsFLogLik)[argNames]) )	#will replace reference to list by contents of list
			# may change parameters on the fly, for permannent changes refer to a list defined in parent and change this one				
			ccl$argsFLogLik = fcl$argsFLogLik
		}
		if( !hasArg("fLogLikScale") && !is.null(fcl$fLogLikScale) ){
			fLogLikScale <- eval.parent(fcl$fLogLikScale)
			ccl$fLogLikScale = fcl$fLogLikScale
		}
	}
}






	




