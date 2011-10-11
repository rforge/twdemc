.calcTemperatedDiffLogDen <- function(
	### calculates the temperated difference of LogDens proposed-accepted per population
	diffLogDenPops		##<< the original differences (comp x pops)
	,TFix				##<< named numeric vector: components with fixed temperature
	,T0c				##<< the temperature per population
){
	diffLogDenPopsT <- diffLogDenPops
	posTFix <- match(names(TFix),rownames(diffLogDenPopsT))
	diffLogDenPopsT[posTFix,,] <- diffLogDenPopsT[posTFix,,,drop=FALSE] / TFix
	for( iPop in 1:dim(diffLogDenPops)[3] ) 
		diffLogDenPopsT[-posTFix,,iPop] <- diffLogDenPopsT[-posTFix,,iPop,drop=FALSE] / T0c[iPop]
	diffLogDenPopsT
}

calcDEMCTempProp <- function(
	### Calculate Temperature of components 
	temp	##<< the maximum temperature
	,diffLogDen		##<< expected difference in LogDensitys proposed-accepted per datastream 
	,rFracMin=1/4	##<< fraction of max DiffDensity  below which temperatue is scaled down to yield  larger importance
){
	#rr <- diffLogDen/max(min(diffLogDen),1e-8)	# diffLogDen per largest misfit (lowest neg diff-LogDensity)
	#rr <- diffLogDen/min(diffLogDen)	# diffLogDen per largest misfit (lowest neg diff-LogDensity)
	rr <- max(diffLogDen/min(diffLogDen), 1e-8)	# ratio of diffLogLik of component per largest abs(diffLogLik), heed that all values are negative, hence min, ratio not below 1e-8
	#Ti <- pmin(temp,1+(temp-1)/rFracMin*rr)
	Ti <- pmin(temp,1+(temp-1)*(rr/rFracMin))
	Ti[rr<=0] <- 1	#give NA or negative values for Ti		
	Ti
	### vector of temperatures corresponding to diffLogDen with maximum corresponding to temp
}

calcDEMCTempDiffLogDen3 <- function(
	### Estimate scalar temperature vector obtain given acceptance rate and scaling temperatures to the same magnitude 
	diffLogDen			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}
	,pTarget=0.18		##<< minimum acceptance rate of component
	,TFix=numeric(0)	##<< named numeric vector: components whose Temperate is to be fixed
	,Tmax=.Machine$double.xmax		##<< maximum temperature to return
	,rFracMin=1/4		##<< fraction of density to which data-streams with low diffLogDen are scaled to by temperatures to
	,doConstrainNeg=FALSE	##<< if given, density of accepted jumps (positive) is constrained to 0
){
	d <- replaceNonFiniteDiffLogDens(diffLogDen, doConstrainNeg=doConstrainNeg)
	l05 <- log(0.5)
	expD <- expDMean <- apply(d,1,median) #rowMeans(d) 
	#apply(d,1,quantile, probs=1-pTarget) #apply(d,1,median)
	iexpD0 <- which(expD>=0) 
	if( 0< length(iexpD0)){
		#will largely understimate temp for those components expD[iexpD0] <- apply(d[iexpD0,,drop=FALSE],1,function(ds){max(ds[ds<0])})
		expD[iexpD0] <- rowMeans(d[iexpD0,,drop=FALSE])
	}
	
	nrowd <- nrow(d)
	Ti <- structure( numeric(nrowd), names=rownames(diffLogDen))
	bo <- (rownames(d) %in% names(TFix))
	Ti[bo] <- TFix
	posNonFix <- (1:nrowd)[!bo]
	
	dVar <- sampleAcceptedFixedTempDiffLogDens(d,TFix=TFix)		#cases of diffLogDen surviving the fixed temp metropolis desicion
	nj <- ncol(dVar)
	
	# objective function to optimize Temperature
	pps <- function(logtemp){  #posNonFix, expD, rFracMin, dTFix, dVar
		temp <- exp(logtemp)
		Ti[posNonFix] <- calcDEMCTempProp( temp, expD[posNonFix], rFracMin )
		dTNonFix <- colSums( dVar[posNonFix,,drop=FALSE]/Ti[posNonFix] )
		#pa <- sum(dTFix+dTNonFix > l05)/nj
		pa <- sum(dTNonFix > l05)/nj
		pa - pTarget
	}
	
	temp2 <- if( ncol(dVar) < 30){
			warning("calcDEMCTempDiffLogDen3: too few accepted cases after fixed temperature rejection, using maximum temperature")
			Tmax
		}else{
			#tmp <- seq(1,1e6,length.out=200) 
			#plot(tmp,pTarget+sapply(tmp,pps))
			if( pps(log(1)) > 0) temp=1  # if acceptance rate at temp 1 is greater than target return temperature 1
			else if( pps(log(Tmax)) < 0) temp=Tmax  	# if acceptance rate at given maximum temperature is already smaller than target return maximum temperature
			else temp <- exp(uniroot( pps, log(c(1, Tmax)), tol=0.01 )$root)
			Ti[posNonFix] <- calcDEMCTempProp( temp, expD[posNonFix], rFracMin )
			TiMean <- Ti	
			
			# in the first round we used the median, do a second round with empirical exptected logDen components: 
			# empirical expected: minimum of actually accepted 
			# based on the Temperatures of the first round
			dT <-  dVar/Ti
			#dTs <- dT[8,]
			pai <- apply(dT,1,function(dTs){ sum(dTs > l05)/nj  })		
			# maximum accepted d per component
			#dTNonFix <- colSums( d[posNonFix,,drop=FALSE]/Ti[posNonFix] )
			boAccepted <- colSums(dT) >= l05
			pAccepted <- sum(boAccepted)/nj
			dAccepted <- dVar[,boAccepted,drop=FALSE]
			expD <- apply(dAccepted,1,min)
			# how many cases would be accepted by using only one components Density ratio:
			# pai <- structure(sapply(seq_along(expD),function(i){sum(d[i,]>=expD[i])/nj} ), names=names(expD))
			if( pps(log(1)) > 0) temp2 = 1 # if acceptance rate at temp 1 is greater than target return temperature 1
			else if( pps(log(Tmax)) < 0) temp2=Tmax  	# if acceptance rate at given maximum temperature is already smaller than target return current temperature
			else temp2 <- exp(uniroot( pps, log(c(1, Tmax)), tol=0.01 )$root)
			temp2
		}
	Ti[posNonFix] <- calcDEMCTempProp( temp2, expD[posNonFix], rFracMin )
	dT <-  dVar/Ti
	pAccepted <- sum(colSums(dT) >= l05)/nj
	attr(Ti,"pAcceptTVar") <- pAccepted
	Ti
	### numeric vector of result components Temperatures
	### yielding acceptance rates closest to pTarget
	### attribute pAcceptTVar giving the calculated acceptance rate for Temperature dependent step
}

calcDEMCTempDiffLogDen3Init <- function(
	### Estimate scalar temperatures to obtain acceptance rate 
	resLogDen			##<< result of \code{\link{twCalcLogDenPar}}: list with components  logDen and resFLogDen
	,...				##<< further arguments to \code{\link{calcDEMCTempDiffLogDenConst}}
	,doConstrainNeg=TRUE	##<< different default value
){
	# replace NA components by sample of non-NA
	# replace non-finite logDens case of smallest finite logDen per component
	for( iComp in 1:ncol(resLogDen$resFLogDen) ){
		boNA <- is.na(resLogDen$resFLogDen[,iComp])
		nNA <- sum(boNA)
		if(0<nNA ){
			resCompNonNA <- resLogDen$resFLogDen[!boNA,iComp]
			resLogDen$resFLogDen[boNA,iComp] <- sample( resCompNonNA, nNA, replace=TRUE)
		}
		boNonFinite <- !is.finite(resLogDen$resFLogDen[,iComp])
		if( 0<sum(boNonFinite)){
			resLogDen$resFLogDen[boNonFinite,iComp] <- min(resLogDen$resFLogDen[!boNonFinite,iComp])
		}
	}
	
	L <- Lp <- resLogDen$resFLogDen 
	#recalculate logDens from replaced components
	rL <- rLQ <- rowSums(resLogDen$resFLogDen)
	
	
	# get the 200 best cases, but at least 5% of the cases
	#rL <- rLFin$logDen
	p <- max(0.05,200/length(rL))
	if( p<1 ){
		q <- quantile(rL, probs=1-p)
		bo <- sapply(q, function(qi) rL >= qi)
		Lp <- L[bo,]
		rLQ <- rL[bo]
	}
	#p1 <- qplot( X2, value, geom="boxplot", data=melt(Lp))+opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p1
	
	#from these 5% sample best as accepted LogDen
	ord <- order(rLQ, decreasing = TRUE)
	La <- Lp[sample( ord[1:ceiling(length(ord)*0.05)], nrow(Lp), replace=TRUE ), ]
	
	#calculate Temperature from density
	d <- Lp-La
	#p2 <- qplot( X2, value, geom="boxplot", data=melt(d))+opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p2
	#T <- calcDEMCTempDiffLogDen(t(d), pTarget=pTarget, TFix=TFix)  # will be scaled in twDEMCInt
	T <- calcDEMCTempDiffLogDen3(t(d), doConstrainNeg=doConstrainNeg, ...)  # will be scaled in twDEMCInt
	T
	### named numeric vector of estimated Temperatures per data stream. 
}

calcDEMCTempDiffLogDen2 <- function(
	### Estimate scalar temperature vector obtain given acceptance rate and scaling temperatures to the same magnitude 
	diffLogDen			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}
	,pTarget=0.18		##<< minimum acceptance rate of component
	,TFix=numeric(0)	##<< named numeric vector: components whose Temperate is to be fixed
	,Tmax=.Machine$double.xmax		##<< maximum temperature to return
	,rFracMin=1/4		##<< fraction of Density to which data-streams with low diff-logDen are scaled to by temperatures to
	,doConstrainNeg=FALSE	##<< if given, density of accepted jumps (positive) is constrained to 0
){
	if( !(is.finite(Tmax) && (Tmax >= 1)) ) stop("need supplying positive TMax")
	#replace non-finite values by lowest finite LogDen value
	#ds <- diffLogDen[1,]
	#diffLogDen <- diffLogDenPops[c("parms","agg_amendm_2","agg_amendm_3","agg_amendm_10"),,1]
	d <- t(apply( diffLogDen,1,function(ds){ds[!is.finite(ds)]<-min(ds[is.finite(ds)]);ds}))
	dimnames(d) <- dimnames(diffLogDen)
	if( doConstrainNeg )
		d[d>0] <- 0
	l05 <- log(0.5)
	expD <- expDMean <- apply(d,1,median) #rowMeans(d) 
	#apply(d,1,quantile, probs=1-pTarget) #apply(d,1,median)
	iexpD0 <- which(expD>=0) 
	if( 0< length(iexpD0)){
		#will largely understimate temp for those components expD[iexpD0] <- apply(d[iexpD0,,drop=FALSE],1,function(ds){max(ds[ds<0])})
		expD[iexpD0] <- rowMeans(d[iexpD0,,drop=FALSE])
	}
	
	nrowd <- nrow(d)
	nj <- ncol(d)
	Ti <- structure( numeric(nrowd), names=rownames(diffLogDen))
	bo <- (rownames(d) %in% names(TFix))
	Ti[bo] <- TFix
	posFix <- (1:nrowd)[bo]
	posNonFix <- (1:nrowd)[!bo]
	
	#x <- d["parms",]
	dTFix <- colSums( d[posFix,,drop=FALSE]/TFix )
	pa <- sum(dTFix>l05)/nj
	if( pa < pTarget ){
		warning("acceptance rate for components with fixed temperature already below target")
		return( Tmax )
	}
	
	pps <- function(logtemp){  #posNonFix, expD, rFracMin, dTFix
		temp <- exp(logtemp)
		Ti[posNonFix] <- calcDEMCTempProp( temp, expD[posNonFix], rFracMin )
		dTNonFix <- colSums( d[posNonFix,,drop=FALSE]/Ti[posNonFix] )
		#apply(d,1,function(ds){sum(ds>l05)/nj})
		pa <- sum(dTFix+dTNonFix > l05)/nj
		pa - pTarget
	}
	
	#tmp <- seq(1,1e6,length.out=200) 
	#plot(tmp,pTarget+sapply(tmp,pps))
	if( pps(log(1)) > 0) temp = 1 # if acceptance rate at temp 1 is greater than target return temperature 1
	else if( pps(log(Tmax)) < 0) temp=Tmax  	# if acceptance rate at given maximum temperature is already smaller than target return current temperature
	else temp <- exp(uniroot( pps, log(c(1, Tmax)), tol=0.01 )$root)
	Ti[posNonFix] <- calcDEMCTempProp( temp, expD[posNonFix], rFracMin )
	TiMean <- Ti	
	
	# in the first round we used the median, do a second round with empirical exptected logDen components: 
	# empirical expected: minimum of actually accepted 
	# based on the Temperatures of the first round
	dT <-  d/Ti
	#dTs <- dT[8,]
	#pai <- apply(dT,1,function(dTs){ sum(dTs > l05)/nj  })
	# maximum accepted d per component
	#dTNonFix <- colSums( d[posNonFix,,drop=FALSE]/Ti[posNonFix] )	
	dAccepted <- d[,colSums(dT) >= l05,drop=FALSE]
	expD <- apply(dAccepted,1,min)
	# how many cases would be accepted by using only one components Density ratio:
	# pai <- structure(sapply(seq_along(expD),function(i){sum(d[i,]>=expD[i])/nj} ), names=names(expD))
	if( pps(log(1)) > 0) temp2 = 1 # if acceptance rate at temp 1 is greater than target return temperature 1
	else if( pps(log(Tmax)) < 0) temp2=Tmax  	# if acceptance rate at given maximum temperature is already smaller than target return current temperature
	else temp2 <- exp(uniroot( pps, log(c(1, Tmax)), tol=0.01 )$root)
	Ti[posNonFix] <- calcDEMCTempProp( temp2, expD[posNonFix], rFracMin )
	
	Ti
	### numeric vector of result components 
	### yielding acceptance rates closest to pTarget
}

.tmp.f.plot.calcDEMCTempDiffLogDen2 <- function(){
	#expD <- seq(-200,0,length.out=101)
	temp=min(expD)/-3
	rr <- expD/min(expD)
	#-------- good Ti <- pmin(temp,1+(temp-1)/rFracMin*rr)
	#-------- exponentially decreasing from Ti to 1: plot( rr, 1+temp-temp^(1-rr) )
	Ti <- pmin(temp,1+(temp-1)/rFracMin*rr)
	#plot(rr,Ti);abline(h=c(1,temp),col="gray");abline(v=rFracMin,col="gray")
	#Ti[1] <- 1
	expDi <- expD/Ti
	plot(expD, expDi)	#for the components below rFracMin the scaled lower temperature increases their weight
	# only for components with a very a diffLogDen close to 0, temperatures do not go below 1 and therefore they have less weight 
	abline(h=min(expDi)*rFracMin, col="gray"); 
	abline(v=min(expD)*rFracMin, col="gray"); 
}

calcDEMCTempDiffLogDen2Init <- function(
	### Estimate scalar temperatures to obtain acceptance rate 
	resLogDen			##<< result of \code{\link{twCalcLogDenPar}}: list with components  logDen and resFLogDen
	,...				##<< further arguments to \code{\link{calcDEMCTempDiffLogDenConst}}
	,doConstrainNeg=TRUE	##<< different default value
){
	# replace NA components by sample of non-NA
	# replace non-finite logDens case of smallest finite logDen per component
	for( iComp in 1:ncol(resLogDen$resFLogDen) ){
		boNA <- is.na(resLogDen$resFLogDen[,iComp])
		nNA <- sum(boNA)
		if(0<nNA ){
			resCompNonNA <- resLogDen$resFLogDen[!boNA,iComp]
			resLogDen$resFLogDen[boNA,iComp] <- sample( resCompNonNA, nNA, replace=TRUE)
		}
		boNonFinite <- !is.finite(resLogDen$resFLogDen[,iComp])
		if( 0<sum(boNonFinite)){
			resLogDen$resFLogDen[boNonFinite,iComp] <- min(resLogDen$resFLogDen[!boNonFinite,iComp])
		}
	}
	
	L <- Lp <- resLogDen$resFLogDen 
	#recalculate logDens
	rL <- rLQ <- rowSums(resLogDen$resFLogDen)
	
	
	# get the 200 best cases, but at least 5% of the cases
	#rL <- rLFin$logDen
	p <- max(0.05,200/length(rL))
	if( p<1 ){
		q <- quantile(rL, probs=1-p)
		bo <- sapply(q, function(qi) rL >= qi)
		Lp <- L[bo,]
		rLQ <- rL[bo]
	}
	#p1 <- qplot( X2, value, geom="boxplot", data=melt(Lp))+opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p1
	
	#from these sample from the 5% best as accepted LogDen
	ord <- order(rLQ, decreasing = TRUE)
	La <- Lp[sample( ord[1:ceiling(length(ord)*0.05)], nrow(Lp), replace=TRUE ), ]
	
	#calculate Temperature from density
	d <- Lp-La
	#p2 <- qplot( X2, value, geom="boxplot", data=melt(d))+opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p2
	#T <- calcDEMCTempDiffLogDen(t(d), pTarget=pTarget, TFix=TFix)  # will be scaled in twDEMCInt
	T <- calcDEMCTempDiffLogDen2(t(d), doConstrainNeg=doConstrainNeg, ...)  # will be scaled in twDEMCInt
	T
	### named numeric vector of estimated Temperatures per data stream. 
}

calcDEMCTempDiffLogDenConst <- function(
	### Estimate scalar temperature to obtain given acceptance rate 
	diffLogDen			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}
	,pTarget=0.18		##<< minimum acceptance rate of component
	,TFix=numeric(0)	##<< named numeric vector: components whose Temperate is to be fixed
	,Tmax=.Machine$double.xmax		##<< maximum temperature to return
){
	#replace non-finite values by lowest finite LogDen value
	#ds <- diffLogDen[1,]
	d <- t(apply( diffLogDen,1,function(ds){ds[!is.finite(ds)]<-min(ds[is.finite(ds)]);ds}))
	l05 <- log(0.5)
	
	nrowd <- nrow(d)
	dnf=d	#non fixed components
	nj <- ncol(d)
	boFix <- TRUE
	if( 0<length(TFix)){
		if( 0 == length(names(TFix))) stop("calcDEMCTempDiffLogDenConst: TFix component must be named")
		TFixPos <- match( names(TFix), rownames(d) )
		if( any(is.na(TFixPos)) ) warning("calcDEMCTempDiffLogDenConst: not all components of TFix in rownames(diffLogDen)")
		dnf<- d[-TFixPos, ,drop=FALSE]
		dfT <- d[TFixPos, ,drop=FALSE]*TFix
		boFix <- apply( dfT,2,function(dfTj){ all(dfTj > l05)})
		#boFix <- rep(TRUE,nj)	#for debugging acceptance rate of others
		#boFix <- sample(c(TRUE,FALSE),nj,replace=TRUE)	#for debugging acceptance rate of others
		
		# check if fixed Temperature components acceptance rate is below pTarget
		paf <- sum(boFix)/nj
		if( paf < pTarget ){
			warning("already fixed Temp components yield acceptance rate below pTarget")
			return(Tmax)
		}
	}
	pps <- function(logtemp){  #dnf, boFix, l05=log(0.5)){	#i: chain
		temp <- exp(logtemp)
		boT <- apply( dnf,2,function(dnfj){all(dnfj/temp > l05)})
		pa <- sum( boFix & boT)/nj
		pa - pTarget 
	}
	#tmp <- seq(1,1e6,length.out=200) 
	#plot(tmp,pTarget+sapply(tmp,pps))
	if( pps(log(1)) > 0) return(1)			# if acceptance rate at temp 1 is greater than target return temperature 1
	if( pps(log(Tmax)) < 0) return(Tmax)  # if acceptance rate at given maximum temperature is already smaller than target return current temperature
	tempOpt <- exp(uniroot( pps, log(c(1, Tmax)), tol=0.01 )$root)
	#exp(uniroot( pps, log(c(1, .Machine$double.xmax)), tol=0.01 )$root)
	# sum(d < log(runif(length(d)))
	tempOpt
	### numeric scalar Temperature between in [1,Tmax] 
	### yielding acceptance rates closest to pTarget
}

calcDEMCTempDiffLogDenConstInit <- function(
	### Estimate scalar temperatures to obtain acceptance rate 
	resLogDen			##<< result of \code{\link{twCalcLogDenPar}}: list with components  logDen and resFLogDen
	,...				##<< further arguments to \code{\link{calcDEMCTempDiffLogDenConst}}
){
	# replace non-finite logDens case of smallest finite logDen per component
	boFin <- is.finite(resLogDen$logDen)
	# replace NA components by sample of non-NA
	for( iComp in 1:ncol(resLogDen$resFLogDen) ){
		resComp <- resLogDen$resFLogDen[,iComp]
		boNA <- is.na(resComp)
		resCompNonNA <- resComp[!boNA]
		nNA <- sum(boNA)
		if(0<length(nNA) )
			resLogDen$resFLogDen[boNA,iComp] <- sample( resCompNonNA, nNA)
	}
	
	boNA <- is.na(resLogDen$l)
	minFiniteLogDen <- min(resLogDen$logDen[boFin])
	minFiniteResFLogDen <- apply( resLogDen$resFLogDen[boFin,], 2, min)
	rLFin <- resLogDen
	rLFin$logDen[!boFin] <- minFiniteLogDen
	for( j in 1:ncol(rLFin$resFLogDen) )
		rLFin$resFLogDen[!boFin,j] <- minFiniteResFLogDen[j]
	
	# get the 200 best cases, but least 30% of the cases
	rL <- rLFin$logDen
	p <- min(max(0.3,200/length(rL)),length(rL))
	q <- quantile(rL, probs=1-p)
	bo <- sapply(q, function(qi) rL > qi)
	Lp <- rLFin$resFLogDen[bo,]
	rLQ <- rLFin$logDen[bo]
	#p1 <- qplot( X2, value, geom="boxplot", data=melt(Lp))+opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p1
	
	#from these sample from the 5% best as accepted LogDen
	ord <- order(rLQ, decreasing = TRUE)
	La <- Lp[sample( ord[1:round(length(ord)*0.05)], nrow(Lp), replace=TRUE ), ]
	
	#calcualte Temperature from density
	d <- Lp-La
	#p2 <- qplot( X2, value, geom="boxplot", data=melt(d))+opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p2
	#T <- calcDEMCTempDiffLogDen(t(d), pTarget=pTarget, TFix=TFix)  # will be scaled in twDEMCInt
	T <- calcDEMCTempDiffLogDenConst(t(d), ...)  # will be scaled in twDEMCInt
	T
	### named numeric vector of estimated Temperatures per data stream. 
}

calcDEMCTempDiffLogDen <- function(
	### Estimate temperatures for different data streams to obtain given acceptance rate 
	diffLogDen			##<< array( streams x steps x  chains) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}
	,pTarget=0.18		##<< overall acceptance rate
	,TFix=numeric(0)	##<< named numeric vector: components whose Temperate is to be fixed 
){
	#replace non-finite values by lowest finite LogDen value
	#ds <- diffLogDen[1,]
	d <- apply( diffLogDen,1,function(ds){ds[!is.finite(ds)]<-min(ds[is.finite(ds)]);ds})
	
	#tmp <- melt(d)
	#p2 <- qplot( X2, value, geom="boxplot", data=tmp)+	opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p2
	
	nrowd <- nrow(d)
	acc <- rep(TRUE, nrowd)
	dq <- d
	if( 0<length(TFix)){
		TFixPos <- match( names(TFix), colnames(d) )
		#TFixPos <- 1:2
		#fpos=1
		tmp1 <- lapply( TFixPos, function(fpos){ di=d[,fpos]/TFix[fpos]; (di > log(runif(nrowd))) })
		for( i in seq_along(tmp1)) acc <- acc & tmp1[[i]]
		# acc now holds rejection by fixed Temperature
		dq <- d[,-TFixPos]
	}
	
	#ps <- 0.8
	ps=0.9
	pps <- function(ps){	#dq,nrowd,acc,pTarget
		qps <- apply( dq, 2, quantile, probs=1-ps, na.rm=TRUE)
		#partial sorting is already efficient qps <- apply( dq, 2, quantileSorted, probs=1-ps, na.rm=TRUE)
		nacc <- sum(sapply(1:nrowd, function(i){ acc[i] & all(dq[i,] > qps) }))
		nacc/nrowd - pTarget
	}
	
	# tmp <- seq(0.2,1,by=0.05)
	# tmp2 <- sapply(tmp, pps)
	# plot( tmp, tmp2 )
	#pps(0.82,0.2)
	ps <- uniroot( pps, c(pTarget, 1), tol=0.001 )$root
	qps <- apply( dq, 2, quantile, probs=1-ps)
	T2q <- pmax(1,qps/log(ps))
	if( 0<length(TFix) ){
		T2 <- structure( numeric(ncol(d)), names=colnames(d) )
		T2[TFixPos] <- TFix
		T2[-TFixPos] <- T2q
	}else 
		T2 <- structure( T2q, names=colnames(d) )
	T2
	### numeric vector of Temperatures (fLogDen_Component x chains)
}

calcDEMCTempDiffLogDenInit <- function(
	### Estimate Temperatures for different data streams to obtain acceptance rate 
	resLogDen			##<< result of \code{\link{twCalcLogDenPar}}: list with components  logDen and resFLogDen
	,...				##<< further arguments to \code{\link{calcDEMCTempDiffLogDen}}
){
	# replace non-finite logDens by very small logDens
	boFin <- is.finite(resLogDen$logDen)
	minFiniteLogDen <- min(resLogDen$logDen[boFin])
	minFiniteResFLogDen <- apply( resLogDen$resFLogDen[boFin,], 2, min)
	rLFin <- resLogDen
	rLFin$logDen[!boFin] <- minFiniteLogDen
	for( j in 1:ncol(rLFin$resFLogDen) )
		rLFin$resFLogDen[!boFin,j] <- minFiniteResFLogDen[j]
	
	# get the 100 best cases, but minium 5% of the cases
	rL <- rLFin$logDen
	p <- min(0.5,100/length(rL),length(rL))
	q <- quantile(rL, probs=1-p)
	bo <- sapply(q, function(qi) rL > qi)
	Lp <- rLFin$resFLogDen[bo,]
	rLQ <- rLFin$logDen[bo]
	#p1 <- qplot( X2, value, geom="boxplot", data=melt(Lp))+opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p1
	
	#from these sample from the 5% best as accepted LogDen
	ord <- order(rLQ, decreasing = TRUE)
	La <- Lp[sample( ord[1:round(length(ord)*0.05)], nrow(Lp), replace=TRUE ), ]
	
	#calcualte Temperature from density
	d <- Lp-La
	#p2 <- qplot( X2, value, geom="boxplot", data=melt(d))+opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p2
	#T <- calcDEMCTempDiffLogDen(t(d), pTarget=pTarget, TFix=TFix)  # will be scaled in twDEMCInt
	T <- calcDEMCTempDiffLogDen(t(d), ...)  # will be scaled in twDEMCInt
	T
	### named numeric vector of estimated Temperatures per data stream. 
}

calcDEMCTemp <- function( 
	### Calculates the temperature for an exponential decrease from \code{T0} to \code{Tend} after \code{nGen} steps. 	
	T0			##<< the initial temperature (before the first step at iGen=0)
	, Tend=1	##<< the temperature at the last step
	, nGen		##<< the number of genrations	
	, iGen=1:nGen ##<< the steps for which to calculate the Temperature	
){
	# calcDEMCTemp
	##seealso<< 
	## \code{\link{twDEMCInt}}
	b = T0
	a = log(Tend/T0)/nGen
	b*exp( a*iGen )
	### vector of Temperatures corresponding to steps iGen
}
attr(calcDEMCTemp,"ex") <- function(){
	plot( 1:100, calcDEMCTemp(T0=100,Tend=5,nGen=100) )	
}

calcDEMCnGenBurnin <- function(
	### Number of steps based on expapolation of an observed temperature decrease to 1 	
	T0 		##<< the initial temperature (before the first step at iGen=0)
	,Ti		##<< the temperatuer at step iStep
	,iStep	##<< the step at which observed Ti
){
	#Ti = T0 e^(a iStep)
	a = log(Ti/T0)/iStep
	-log(T0)/a
}

calcDEMCTempGlobal1 <- function(
	### Calculating global temperature after the next batch for one population.
	resPop		##<< twDEMC result (subChains of one population)
	,diffLogDen	##<< numeric vector: Lp-La of the previous proposals
	,TLp		##<< numeric scalar: max Temperature suggested by optimizing Lp  (from \code{\link{calcDEMCTempDiffLogDen3}}			
	,pAcceptTVar ##<< numeric scalar: Acceptance rate of temperatue dependent step (from \code{\link{calcDEMCTempDiffLogDen3}}
	,iRun=getNGen(resPop)	##<< current generation: may be passed for efficiency
	,nGenBurnin ##<< integer scalar: the number of Generations in burnin
	,nRun		##<< integer scalar: the number of generations in next batch
	,minPCompAcceptTempDecr=0.16
){
	##details<< 
	## This version either enforces Temp-Decrease complying to exponential decrease to 1 at nGenBurnin
	## or stays at temperature and prolonges nGenBurnin
	T0 <- resPop$temp[nrow(resPop$temp),1]
	if( pAcceptTVar < minPCompAcceptTempDecr){
		TGlobal <- T0
		nGenBurnin <- nGenBurnin+nRun
	}else{
		TGlobal <- calcDEMCTemp( T0, 1, nGenBurnin-iRun, nRun)
	} 
	list(TGlobal=TGlobal,nGenBurnin=nGenBurnin)	
	### list with components \itemize{
	### \item{TGlobal: numeric scalar: the global Temperature}
	### \item{nGenBurnin: recalculated burnin period}
	### }
}

calcDEMCTempGlobal2a <- function(
	### Calculating global temperature after the next batch for one population.
	resPop		##<< twDEMC result (subChains of one population)
	,diffLogDen	##<< numeric vector: Lp-La of the previous proposals
	,TLp		##<< numeric scalar: max Temperature suggested by optimizing Lp  (from \code{\link{calcDEMCTempDiffLogDen3}}			
	,pAcceptTVar ##<< numeric scalar: Acceptance rate of temperatue dependent step (from \code{\link{calcDEMCTempDiffLogDen3}}
	,iRun=getNGen(resPop)	##<< current generation: may be passed for efficiency
	,nGenBurnin ##<< integer scalar: the number of Generations in burnin
	,nRun		##<< integer scalar: the number of generations in next batch
	,rHat0=1.08	##<< rHat value for which to not change the burnin period
){
	rHat0=max(rHat0,1.1)		#not smaller than 1.1 otherwise too strong  
	nR <- nrow(resPop$temp)
	T0 <- resPop$temp[nR,1]
	# gelman diagnostics for the last part of Temperature decrease from 120% of current T
	i130 <- twBinOptimize(resPop$temp[,1], 1.2*T0)$where
	i130b <- max(1,min(nR-30,i130)) # go at least 100 rows back (but not over beginning)
	if( resPop$temp[i130b,1] < 1.2*T0 ){
		nGenBurninNew <- nGenBurnin		# if Temperature did not change do not adjust interval
		TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
	}else{
		res130 <- thin(resPop, start=i130b)
		tmp <- checkConvergenceGelman(res130,burninFrac=0)
		r2 <- max(attr(tmp,"rHat"))
cat("r2=",r2,"\n",sep="")		
		# map gelman diag to a change in nGenBurnin (multiply the period until nGenBurnin by a factor)
		#r2 <- seq(0.9,2,by=0.02)
		burninFac <- pmax(-2/3, log(1/3)/(rHat0-1.0) * (rHat0-r2))
		#plot(burninFac~r2); abline(h=1); abline(v=rHat0)
		if( burninFac > 0){
			# do not decrease Temperature in the next run
			TGlobal <- T0
			# increase by factor + next batch
			nGenBurninNew <- round(nGenBurnin+(nGenBurnin-iRun)*burninFac)+nRun
		}else{
			# decrease burnin period
			nGenBurninNew <- round(nGenBurnin+(nGenBurnin-iRun)*burninFac)
			TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
		}
	}
	#i <- iRun:(nGenBurninNew*1.2)
	#rT <- T0^(1/(nGenBurninNew-iRun))
	#plot( rT^(nGenBurninNew-i) ~ i)
	#TGlobal <- rT^(nGenBurninNew-(iRun+nRun))
	list(TGlobal=TGlobal,nGenBurnin=nGenBurninNew)	
	### list with components \itemize{
	### \item{TGlobal: numeric scalar: the global Temperature}
	### \item{nGenBurnin: recalculated burnin period}
	### }
}

calcDEMCTempGlobal2b <- function(
	### Calculating global temperature after the next batch for one population.
	resPop		##<< twDEMC result (subChains of one population)
	,diffLogDen	##<< numeric vector: Lp-La of the previous proposals
	,TLp		##<< numeric scalar: max Temperature suggested by optimizing Lp  (from \code{\link{calcDEMCTempDiffLogDen3}}			
	,pAcceptTVar ##<< numeric scalar: Acceptance rate of temperatue dependent step (from \code{\link{calcDEMCTempDiffLogDen3}}
	,iRun=getNGen(resPop)	##<< current generation: may be passed for efficiency
	,nGenBurnin ##<< integer scalar: the number of Generations in burnin
	,nRun		##<< integer scalar: the number of generations in next batch
	,rHat0=1.06	##<< rHat value for which to not change the burnin period
){
	rHat0=max(rHat0,1.1)		#not smaller than 1.1 otherwise too strong
	temp <- resPop$temp[,1]
	nR <- length(temp)
	T0 <- temp[nR]
	acceptRowsFac <- 1/(resPop$thin*resPop$pAccept[nR,1])
	i <- max(1,round(nR- 20*acceptRowsFac))	# row 20 independent steps back
	t20r<-temp[i]/T0
	#if( iRun <= nRun ){
		# assume that nRun was also the previous batch
		# after first batch do not adjust temperature
		# if Temperature did not change do not adjust interval and decrease temperature
	#	cat("  first batch: no adjustment of burnin length\n")		
	#	nGenBurninNew <- nGenBurnin		
	#	TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
	#}else{
		minRTback20ToTCurr <- 1.25
		if( t20r > minRTback20ToTCurr){
			# 20 accepted rows back, temperate was more than 125% of current temperate: too fast
			# keep current temperature and extend burnin by generations in next batch
			cat("  Tback20/TCurr=",round(t20r*100),"%, pAccept=",resPop$pAccept[nR,1],"\n",sep="")		
			TGlobal <- T0
			nGenBurninNew <- nGenBurnin*1.1+nRun
		}else{
			#i130 <- twBinOptimize(temp[1:(iRun-nRun)], 1.2*T0)$where	# here do not regard the last batch
			#if( is.na(i130) || (temp[i130] < 1.1*T0) ){
		#		# if Temperature did not change do not adjust interval and decrease temperature
		#		nGenBurninNew <- nGenBurnin		
	    #			TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
	    #	}else {
			i130 <- twBinOptimize(temp[temp!=T0], 1.2*T0)$where	
			if(is.na(i130)){
				cat("  no upper temperature found: no adjustment of burnin length\n")		
				nGenBurninNew <- nGenBurnin		
				TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
			}else{
				tAccRows120 <- (iRun-i130)/acceptRowsFac 
				if( tAccRows120 < 15 ){
					#temperature decrease from 120% to 100% in less than 15 accepted rows: too fast
					# keep current temperature and extend burnin by generations in next batch
					cat("  accepted rows from a 20% Temperature decrease=",tAccRows120,", pAccept=",resPop$pAccept[nR,1],"\n",sep="")		
					TGlobal <- T0
					nGenBurninNew <- nGenBurnin*1.1+nRun
				}else{
					i130b <- max(1,i130) # go at least 15 accepted steps back (but not over beginning)
					#if( temp[i130b] < 1.2*T0 ){
					#	cat("  no temperature change within 15 accepted steps: no adjustment of burnin length\n")		
					#	nGenBurninNew <- nGenBurnin		# if Temperature did not change do not adjust interval
					#	TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
					#}else{
						# Gelman diag on properly thinned population
						#dump.frames(file.path("tmp","tempDecGelman"),TRUE)
						#stop("dump to tmp/tempDecGelman.rda")
						newThinOdd <- 1/resPop$pAccept[nR,1]
						newThin <- max(1,(newThinOdd%/%resPop$thin))*resPop$thin	# make it multiple of current thin
						res130 <- thin(resPop, start=i130b, newThin=newThin)
						tmp <- checkConvergenceGelman(res130,burninFrac=0)
						r2 <- max(attr(tmp,"rHat"))
						cat("  Gelman diag: r2=",r2,"\n",sep="")		
						# map gelman diag to a change in nGenBurnin (multiply the period until nGenBurnin by a factor)
						#r2 <- seq(0.9,2,by=0.02)
						burninFac <- pmax(-2/3, log(1/3)/(rHat0-1.0) * (rHat0-r2))
						#plot(burninFac~r2); abline(h=1); abline(v=rHat0)
						if( burninFac > 0){
							# do not decrease Temperature in the next run
							TGlobal <- T0
							# increase by factor + next batch
							nGenBurninNew <- round(nGenBurnin+(nGenBurnin-iRun)*burninFac)+nRun
						}else{
							# decrease burnin period
							nGenBurninNew <- round(nGenBurnin+(nGenBurnin-iRun)*burninFac)
							TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
						}
					} # Gelman Diag
				#} #120% to 100% in less than 15 accepted rows
			} # no upper temperature found
		} #20 accepted rows back
	#} # first batch
	list(TGlobal=TGlobal,nGenBurnin=nGenBurninNew)	
	### list with components \itemize{
	### \item{TGlobal: numeric scalar: the global Temperature}
	### \item{nGenBurnin: recalculated burnin period}
	### }
}

.tmp.f <- function(){
	load("../asom/tmp/tempDecGelman.rda")
	debugger("tmp/tempDecGelman")
	load("../tmp/tempDecGelman.rda.rda")
	debugger(get("tmp/tempDecGelman.rda"))
}

.tmp.f <- function(){
	# T(rT,i) = rT^((nGenBurnin-nRun)-i)
	# T(rT,0) = T0
	rT <- T0^(1/(nGenBurnin+0-iRun))
	rT2 <- T0^(1/(nGenBurnin+(nGenBurnin-nRun)*burninFac-nRun))
	rT3 <- T0^(1/(nGenBurnin+(nGenBurnin-nRun)*1-nRun))
	rT4 <- T0^(1/(nGenBurnin+(nGenBurnin-nRun)*-0.5-nRun))
	i<-nRun:nGenBurnin
	plot(rT^(nGenBurnin-i)~i); 
	lines( rT2^(nGenBurnin+(nGenBurnin-nRun)*burninFac-i)~i,col="blue" ) 
	lines( rT3^(nGenBurnin+(nGenBurnin-nRun)*1-i)~i,col="red" )
	lines( rT4^(nGenBurnin+(nGenBurnin-nRun)*-0.5-i)~i,col="maroon" )
	
}

#mtrace(calcDEMCTempGlobal2)
calcDEMCTempGlobal2 <- function(
	### Calculating global temperature and adjusted burnin period after the next batch for one population.
	resPop		##<< twDEMC result (subChains of one population)
	,diffLogDen	##<< numeric vector: Lp-La of the previous proposals
	,TLp		##<< numeric scalar: max Temperature suggested by optimizing Lp  (from \code{\link{calcDEMCTempDiffLogDen3}}			
	,pAcceptTVar ##<< numeric scalar: Acceptance rate of temperatue dependent step (from \code{\link{calcDEMCTempDiffLogDen3}}
	,iRun=getNGen(resPop)	##<< current generation: may be passed for efficiency
	,nGenBurnin ##<< integer scalar: the number of Generations in burnin
	,nRun		##<< integer scalar: the number of generations in next batch
	,rHat0=1.2	##<< rHat value for which to not change the burnin period
){
	rHat0=max(rHat0,1.001)		#not smaller than 1.001 in order to avoid division by zero
	temp <- resPop$temp[,1]
	nR <- length(temp)
	T0 <- temp[nR]
	acceptRowsFac <- 1/(resPop$thin*resPop$pAccept[nR,1])  # every x rows can be regarded as independent
	i <- max(1,round(nR- 20*acceptRowsFac))	# row 20 independent steps back
	t20r<-temp[i]/T0
	if( t20r > 1.2){
		# 20 accepted rows back, temperate was more than 120% of current temperate: too fast
		# keep current temperature and extend burnin by generations in next batch
		cat("  Tback20/TCurr=",round(t20r*100),"%, pAccept=",resPop$pAccept[nR,1]," keep T0 and extend burnin\n",sep="")		
		TGlobal <- T0
		nGenBurninNew <- nGenBurnin*1.1+nRun
	}else{
		temp2 <- temp[temp!=T0]
		i130 <- if(length(temp2)>0) twBinOptimize(temp2, 1.2*T0)$where else NA
		if(is.na(i130) || (temp[i130]<1.1*T0) ){
			cat("  no significant temperature decrease yet: no adjustment of burnin length\n")		
			nGenBurninNew <- nGenBurnin		
			TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
		}else{
			newThinOdd <- 1/resPop$pAccept[nR,1]
			newThin <- max(1,(newThinOdd%/%resPop$thin))*resPop$thin	# make it multiple of current thin
			res130 <- thin(resPop, start=max(i130, iRun-2*nRun), newThin=newThin)
			# if the chains of one populatin cover the minimum at given temperature, temperature may decrease
			# hence check, if they converged to limiting distribution for given temperature
			if( checkConvergenceTrend(res130) < 0.05){
				# for a significant trend in rLogDen of resPop, stay at given temperature
				cat("  trend in logDensity: stay at given temperature one more batch\n")		
				nGenBurninNew <- nGenBurnin+nRun
				TGlobal <- T0
			}else{
				# base diagnostics on end of the chain with a temperatue decrease of 10%
				#dump.frames(file.path("tmp","tempDecGelman"),TRUE)
				#stop("dump to tmp/tempDecGelman.rda")
				d <- dim(res130$parms)[2:3]
				xGrid <- rep(1:d[1],d[2])
				rHat2 <- apply(res130$parms,1,function(parmsi){
					# before Gelman diag, first remove trend common to all chains 
					# that will otherwise dominate both variances
					# with decreasing temperature a trend is very probable, despite the Gelman diag measures the mixing of the chains
					parmsid <- .detrendMatrix(parmsi,xGrid=xGrid,df=3)
					rHat2 <- twDEMC:::.calcRhat2( parmsid, n=d[1], m=d[2] )
				})
				rHatMax <- sqrt(max(1,rHat2))
				cat("  Gelman diag: max(rHat)=",rHatMax,": ",ifelse(rHatMax <rHat0,"decreasing","increasing")," burnin\n",sep="")		
				# map gelman diag to a change in nGenBurnin (multiply the period until nGenBurnin by a factor)
				#rHatMax <- seq(1.0,2,by=0.02)
				#burninFac <- pmax(-1/2, log(1/3)/(rHat0-1.0) * (rHat0-rHatMax))	# decrease to strong
				burninFac <- pmax(-1/2, log(1/2)/(rHat0-1.0) * (rHat0-rHatMax))
				#plot(burninFac~rHatMax); abline(h=0); abline(v=rHat0)
				#rHat0 <- seq(1.00,1.4,by=0.01); 
				#plot( burninFac ~ rHat0 ); abline(h=c(log(1/3),0,1));
				if( burninFac > 0){
					# do not decrease Temperature in the next run
					TGlobal <- T0
					# increase by factor + next batch
					nGenBurninNew <- round(nGenBurnin+(nGenBurnin-iRun)*burninFac)+nRun
					#plot( nGenBurninNew ~ burninFac ); abline(h=nGenBurnin); abline(v=0)
					#plot( nGenBurninNew ~ rHatMax ); abline(v=rHat0); abline(h=nGenBurnin)
				}else{
					#decrease Temperature and shorten burnin
					nGenBurninNew <- (nGenBurnin+(nGenBurnin-iRun)*burninFac)
					TGlobal <- T0^(1-nRun/(nGenBurninNew-iRun))
				} 
			}# end Gelman Diag 
		} # end no significnat temperature found
	} #20 end accepted rows back
	# make sure that Temperature decrease is not gerater than 100/120 within next 20  accepted steps (a)
	a <- 20/(0.8*resPop$pAccept[nR,1])	# calculate with slightly lower acceptance rate to get a robust estimate
	Ta = 1.2^(1/a)				# T and burnin so that Ta^(ba-(i+a) = 100/120 T0
	ba = round(log(T0)/log(Ta)+iRun)
	TaGlobal <- Ta^(ba-(iRun+nRun))
	TGlobal <- max(TaGlobal,TGlobal)
	nGenBurninNew <- max(ba,round(nGenBurninNew))

	# calculate T
	list(TGlobal=TGlobal,nGenBurnin=nGenBurninNew)	
	### list with components \itemize{
	### \item{TGlobal: numeric scalar: the global Temperature}
	### \item{nGenBurnin: recalculated burnin period}
	### }
}

.tmp.f <- function(){
	# exploring Temperature decrease of 100/120 after a generations
	b=100
	iRun=20
	a=20
	i= 1:b
	T0=10
	T = T0^(1/(b-iRun))
	plot( T^(b-i) ~ i);	abline(h=T0); abline(v=iRun)
	
	Ta = 1.2^(1/a)
	ba = log(T0)/log(Ta)+iRun
	lines( Ta^(ba-i) ~ i,col="blue");
	T0a = Ta^(ba-(iRun+a))
	T0/T0a
}

.tmp.f <- function(){
	# exploring smooting splines on vairable
	pa <- res130$parms["a",,]
	xGrid <- 1:nrow(pa)
	(spl <- smooth.spline(rep(xGrid,ncol(pa)),as.vector(pa),df=3 ))
	matplot(pa,type="p")
	lines(spl$y~xGrid)
	
	pa2 <- pa - spl$y
}

.detrendMatrix <- function(
	### removes a common trend (smoothed data) across all columns
	pa					##<< numeric matrix to be detrended
	,df=3 				##<< degrees of freedeom for the trend, see \code{\link{smooth.spline}}
	,xGrid=rep(1:nrow(pa),ncol(pa))	##<< for several similar matrices, may be passed for performance reasons
){
	trend <- smooth.spline(xGrid,as.vector(pa),df=df )$y
	pa - trend
	### numeric matrix 
}


