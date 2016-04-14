
ofMultiSumArgsList <- function( 
	### sum of objective function returning vector results (multiple data streams)
	popt,					##<< the parameter vector at which ofMulti is evaulated 
	argsList, 				##<< argument list for the objective function, first component must be the objective function itself
	returnComponents=FALSE,	##<< see return value  
	#prevLogDenParms=NULL,  	##<< previous logDenlihood argument appended to argsList
	... 					##<< further arguments passed to objective function
){
	# ofMultiSumArgsList
	
	##details<< 
	## If argsList is a name, it is first evaulated in the parent frame.
	## This supports the usage of sfExport of the argsList instead of passing much data in sfApply via \dots
	
	# #details<< 
	# # prevLogDenParms is supplied to ofMulit in order to support 
	# # a two-step Metropolis to avoid unneccesary model evaluations (-1/2 misfit)
	
	if( is.name(argsList) ) argsList = eval.parent(argsList)
	body <- expression({
			#cargs = c(list(argsList[[1]]), list(popt), argsList[-1])
			cargs = c(list(popt), argsList[-1], list(...) )
			# if( !is.null(prevLogDenParms) ) cargs = c( cargs, list(prevLogDenParms=prevLogDenParms) )
			#(comp <- eval(as.call(cargs)))
			(comp <- do.call(argsList[[1]], cargs))
		})
	##details<<  
	## The evaluated argsList is checked for entry dumpFileBaseName. If it exists
	## a dump is created in this file when an error occurs and execution stops with an error
	if( !is.null(argsList$dumpFileBaseName)){
		dumpFileBaseName <- argsList$dumpFileBaseName; argsList$dumpFileBaseName <- NULL	#do not provide argument to fLogDen
		comp <- try( eval(body) )
		if( inherits(comp, "try-error") ){
			dump.frames(dumpFileBaseName,TRUE)
			stop(comp)
		}
	}else{
		comp <- eval(body)
	}
	#sumcomp = if(is.null(prevLogDenParms)) sum(comp) else sum( comp[ names(comp)!="parms"] )	#assume parms is handled in first Metropolis-Step inside fLogDen
	sumcomp =  sum(comp)	#make sure that with a two-step Metropolis, parms is not returned as one components 
	if( returnComponents ) list(sum=sumcomp, comp=comp ) else sumcomp
	### scalar sum of result of objective function. 
	### If returnComponents is true, value is not a scalar but a list 
	### with the following compoents \describe{
	### \item{sum}{ the scalar sum as before }
	### \item{sumcomp}{ the vector result of the objective function}}
}
#mtrace(ofMultiSumArgsList)
#load("testdump.rda")
#debugger(testdump)



ofMultiSum <- function( 
	### weighted sum of objective function returning vector results (multiple data streams)
	popt,				##<< parameters for which to calculate objective function 
	ofMulti=NULL, 		##<< objective function
	lowerBoundPopt=NULL, ##<< additionally returns a large penalty value if parameters are outside lower or upper bound
	upperBoundPopt=NULL, 
	outsideBoundVal=Inf, 
	multiWeights=NULL,	##<< weights: column names must correspond to the column of the return of ofMulti 
	...					##<< further argumetns passed to ofMulti
){
	# ofMultiSum
	r <- 0	
	if( !is.null(lowerBoundPopt) )
		if( any(popt-lowerBoundPopt < 0)) r <- outsideBoundVal
	if( !is.null(ub) )
		if( any(ub-popt < 0)) r <- outsideBoundVal
	if( r == 0){
		rcc <- ofMulti(popt, ...)
		#rcc = try( ofMulti(popt, ...), TRUE )
		if( inherits(rcc,"try-error")){
			r=Inf
		}else{
			if( is.null(multiWeights) ){
				r =sum( rcc )
			}else{
				r =sum( multiWeights[names(rcc)]*rcc )
			}
		}
	}
	r
}
#mtrace(ofMultiSum)
#sfExport("ofMultiSum")

ofMultiSumSimple <- function( 
	### sum of objective function returning vector results (multiple data streams)
	popt, 
	ofMulti, 
	...
){
	sum(ofMulti(popt, ...))
}
#mtrace(ofMultiSumSimple)

