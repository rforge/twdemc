# XX not adapted to block sampling yet

.tmp.f <- function(){
getAcceptedPos.twDEMCProps <- function(
	### Calculate the accepted state corresponding to the proposal.
	Y	##<< matrix of proposals with row "accepted" and first step (column) initial state, and third dimension chains.
	,acceptedRowName = "accepted"
){
	#j=1
	# set the accepted state of the first entry to TRUE so that a start is given 
	Y[acceptedRowName,1,] <- 1
	tmp.f <- function(j){		#j is the chain, i.e. column
		acc <- which( Y[acceptedRowName,,j,drop=TRUE] != 0)
		accl <- diff(c(acc, ncol(Y)+1))
		rep(acc,accl)
	}	
	acceptedPos <- sapply(1:dim(Y)[3], tmp.f) 
	#acceptPos for an accepted proposal is the proposal step itselv, we need the previous one
	if( is.matrix(acceptedPos) )
		prevAcceptedPos <- acceptedPos[c(1,1:(nrow(acceptedPos)-1)),]
	else
		# sapply reduces the case for only 1 step to a vector
		prevAcceptedPos <- matrix(acceptedPos,nrow=1)
	prevAcceptedPos[1,] <- NA
	prevAcceptedPos
	### matrix of steps (steps x chains) whose state was the accepted state when evaluating proposal (row of Y)
	### the first step is NA because there is no previous accepted state
}

getDiffLogDen.twDEMCProps <- function(
	### Extract the Differences in LogDensity between accepted states and proposals.
	Y					##<< matrix of proposals with row "accepted" and first step (column) initial state, rows: results components of fLogDen and third dimension chains.
	,resCols			##<< the rows of Y with result components, either names or positions
	,nLastSteps = 128	##<< number of last steps of Y for which to extract diffs
	,temp=1				##<< numeric matrix (resComp x pop): the temperature applied to the difference of LogDens
	,...				##<< further arguments passed to \code{\link{getAcceptedPos.twDEMCProps}}
){
	#work with positions
	if( is.character(resCols) ) resCols <- match(resCols,rownames(Y))
	# make sure that temp has all resCols components 
	if( 1 < length(temp)) temp <- temp[rownames(Y)[resCols],] else temp <- matrix(temp,ncol=dim(Y)[3],nrow=1)
	nPopsChain <- dim(Y)[3] %/% ncol(temp)
	# constrain to the nLastSteps+1 rows
	YL <- if( nLastSteps+1 >= ncol(Y)) Y else Y[,ncol(Y)+1-((nLastSteps+1):1),,drop=FALSE]
	#test 1d case: YL <- Y[,min(ncol(Y),nLastSteps+1),,drop=FALSE]
	#mtrace(getAcceptedPos.twDEMCProps)
	acceptedPos <- getAcceptedPos.twDEMCProps(YL,...) #index of parameter vector that correponds to the currently accepted for vector at given row
	diffLogDen <- abind(lapply( 1:dim(YL)[3], function(j){ 	adrop((YL[resCols,-1,j,drop=FALSE] - YL[resCols,acceptedPos[-1,j],j,drop=FALSE])/temp[,(j-1)%/%nPopsChain+1],3)}),rev.along=0)
	diffLogDen
	### numeric array ( component x nLastSteps x chain ) of Lp-La, the L
}


replaceNonFiniteDiffLogDens <- function(
	### For each component replace NAs by sample of others and non-Finite values by minimum of others  
	diffLogDen			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}
	,doConstrainNeg=FALSE	##<< if given, density of accepted jumps (positive) is constrained to 0
){
	d <- t(apply( diffLogDen,1,function(ds){
				boNA <- is.na(ds)
				nNA <- sum(boNA)
				if(0<nNA ){
					dsNonNA <- ds[!boNA]
					ds[boNA] <- sample( dsNonNA, nNA, replace=TRUE)
				}
				boNonFinite <- !is.finite(ds)
				if( 0<sum(boNonFinite)){
					ds[boNonFinite] <- min(ds[!boNonFinite])
				}
				ds
			}))
	dimnames(d) <- dimnames(diffLogDen)
	if( doConstrainNeg )
		d[d>0] <- 0
	d
}

sampleAcceptedFixedTempDiffLogDens <- function(
	### Calculate diffLogDen for fixed temperature components and return subset of accepted ones 
	diffLogDen			##<< array( streams x steps) Lp-La see \code{\link{getDiffLogDen.twDEMCProps}}
	,TFix=numeric(0)
){
	dTFix <- colSums( diffLogDen[names(TFix),,drop=FALSE]/TFix )
	acceptedFixed <- dTFix>log(runif(ncol(diffLogDen)))
	#pa <- sum(acceptedFixed)/nj
	diffLogDen[,acceptedFixed,drop=FALSE]
}

} #.tmp.f

