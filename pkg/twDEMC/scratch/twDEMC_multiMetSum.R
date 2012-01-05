
# failed to base temperature on percentiles
# better use the version with "all", 
calcDEMCTempDiffLogDenConst <- function(
	### Estimate scalar Temperature to obtain given acceptance rate 
	diffLogDen			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCPops}}
	,pTarget=0.2		##<< overall acceptance rate
	,TFix=numeric(0)	##<< named numeric vector: components whose Temperate is to be fixed
	,TCurr=.Machine$double.xmax		##<< current Temperature 
){
	#replace non-finite values by lowest finite LogDen value
	#ds <- diffLogDen[1,]
	d <- t(apply( diffLogDen,1,function(ds){ds[!is.finite(ds)]<-min(ds[is.finite(ds)]);ds}))
	l05 <- log(0.5)
	
	nrowd <- nrow(d)
	dF= numeric(ncol(d)) # sum fixed components*T (steps x chains)
	dnf=d	#non fixed components
	if( 0<length(TFix)){
		if( 0 == length(names(TFix))) stop("calcDEMCTempDiffLogDenConst: TFix component must be named")
		TFixPos <- match( names(TFix), rownames(d) )
		if( any(is.na(TFixPos)) ) warning("calcDEMCTempDiffLogDenConst: not all components of TFix in rownames(diffLogDen)")
		dF <- colSums( d[TFixPos, ,drop=FALSE]*TFix)	
		dnf<- d[-TFixPos, ,drop=FALSE]
		
		# check if fixed Temperature components acceptance rate is below pTarget
		#pa <- 1-ecdf(dF)(l05) #percentil of log(0.5)
		f1 <- ecdf(dF)
		pa <- 1-approx(knots(f1),f1(knots(f1)),xout=l05)$y	#for linear interpolation in between knots
		if( pa < pTarget){
			warning("acceptance rate based on fixed Temp components is already below target level")
			return(TCurr)
		}
	}
	dV <- colSums(dnf)	#sum across non-fixed components (steps x chains)
	
	dF <- sort(dF)
	dV <- sort(dV)
	#temp <- 50 
	pps <- function(temp){  #,dF,dV,pTarget,l05=log(0.5)){	#i: chain
		dT <- dF + dV/temp
		tmpi <- max(which(dT<l05))
		#f1 <- ecdf(dt)
		#pa <- 1-approx(knots(f1),f1(knots(f1)),xout=l05)$y	#for linear interpolation in between knots
		pa <- 1-ecdf(dT)(0.5)	#percentil of log(0.5)
		pa - pTarget 
	}
	
	#tmp <- 100:1
	#plot(tmp,sapply(tmp,pps))
	if( pps(1) > 0) return(1)
	tempOpt <- uniroot( pps, c(1, TCurr), tol=0.001 )$root
	# sum(d < log(runif(length(d))) 
	### numeric scalar Temperature
}

