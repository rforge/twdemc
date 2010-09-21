sfApplyMatrixLB <- function(
	### load balanced parallel application of FUN to subMatrices X1 and X2   
	X,			##<< matrix with rows corresponding to cases supplied to first argument of FUN 
	MARGIN=1,	##<< 1 indicates rows, 2 indicates columns
	FUN,		##<< function to be applied on each row
	...,		##<< further arguments passed to FUN
	debugSequential = FALSE 	##<< use apply 
){
	if( !is.matrix(X) )
		stop("X must be matrix.")
	if( !(MARGIN %in% 1:2))
		stop("MARGIN must be 1 or 2.")
	if( nrow(X) == 0)
		return( list() )
	##details<< 
	## if debugSequentail is TRUE, simple apply is used 
	res <- if( debugSequential ){
			apply(X, MARGIN, FUN, ... )
		}else{
			##details<< if debugSequential is FALSE (default), then  
			## internally a variable with random name is created for X and exported to all slaves.
			## An intermediate function is called on each cluster node provided with the index and the names of the exported variables.
			## The internal function evaluates X again and executes FUN[X1[index,],...) or FUN[X1[,index],...) respectively 
			XName <- paste("X",sample.int(32000,1),sep="_")
			XNameName <- as.name(XName)
			assign(XName,X)
			sfExport(list=as.list(XName))
			applyFunc <- if(MARGIN==1) .twSfApplyMatrixRowsLB_internal else .twSfApplyMatrixColsLB_internal
			nCases <- if(MARGIN==1) nrow(X) else ncol(X)
			caseNames <- dimnames(X)[MARGIN]
			resl <- sfClusterApplyLB( 1:nCases, applyFunc, XNameName, FUN, ... )
			sfClusterCall( function(XName){rm(list=c(XName),pos=1)}, XName )
			#does not work??: sfRemove(list=XNames)
			rm(list=XName)
			res <- .simplifyLBResult(resl,caseNames)
		}#!debugSeqential
	### a list with result of FUN for each row of X1 and X2
	
	##seealso<<   
	## \code{\link{twDEMCInt}}
}
# mtrace(sfApplyMatrixLB)
# twUtestF("applyLB","test.sfApplyMatrixLB")
# twUtestF("applyLB")

#.simplifyLBResult not in parallel.R 

.twSfApplyMatrixRowsLB_internal <- function(
	### internal function called on each cluster node from \code{\link{sfApplyMatrixLB}}
	index,XNameName,FUN,...
){
	X <- eval.parent(XNameName)
	FUN(X[index,],...)
}

.twSfApplyMatrixColsLB_internal <- function(
	### internal function called on each cluster node from \code{\link{sfApplyMatrixLB}}
	index,XNameName,FUN,...
){
	X <- eval.parent(XNameName)
	FUN(X[,index],...)
}



twSfApply2MatrixLB <- function(
	### load balanced parallel application of FUN to subMatrices X1 and X2   
	X1,		##<< matrix with rows corresponding to cases supplied to first argument of FUN 
	X2,		##<< matrix with the same number of rows as X1, cases supplied to second argument of FUN 
	FUN,	##<< function to be applied on each row
	...,
	debugSequential = FALSE,
	simplifyLBResult= TRUE
){
	if( !is.matrix(X1) | !is.matrix(X2))
		stop("X1 and X2 must be matrices.")
	if( !identical(nrow(X1),nrow(X2)) )
		stop("Number of rows in X1 and X2 must be the same.")
	if( nrow(X1) == 0)
		return( list() )
	##details<< 
	## if debugSequentail is TRUE, simple lapply is used 
	resl <- if( debugSequential ){
				lapply(1:nrow(X1), function(index,FUN,...){FUN(X1[index,],X2[index,],...)} ,FUN,... )
		}else{
			##details<< if debugSequential is FALSE (default), then  
			## internally variables with random names are created for X1 and X2 and exported to all slaves.
			## An intermediate function is called on each cluster node provided with the index and the names of the exported variables.
			## The internal function evaluates X1 and X2 again and executes FUN[X1[index,],X2[index,],...)
			XNames <- paste("X",1:2,sample.int(32000,1),sep="_")
			XNamesNames <- lapply(XNames,as.name)
			assign(XNames[1],X1)
			assign(XNames[2],X2)
			sfExport(list=as.list(XNames))
			#sfClusterEval(ls())
			#lapply(1:nrow(X1), .twSfApply2MatrixLB_internal, XNamesNames, paste, sep="_" )
			#sfClusterApplyLB(1:nrow(X1), .twSfApply2MatrixLB_internal, XNamesNames, paste, sep="_" )
			res <- sfClusterApplyLB( 1:nrow(X1), .twSfApply2MatrixLB_internal, XNamesNames, FUN, ... )
			sfClusterCall( function(XNames){rm(list=c(XNames),pos=1)}, XNames )
			#does not work?? sfRemove(list=XNames)
			rm(list=XNames)
			res
		}
	if( simplifyLBResult ) .simplifyLBResult(resl, dimnames(X1)[1]) else resl
	### result of \code{\link{.simplifyLBResult}} applied to a list with result of FUN for each row of X1 and X2 
	
	##seealso<<   
	## \code{\link{twDEMCInt}}
}
# twUtestF("applyLB","test.twSfApply2MatrixLB")

.twSfApply2MatrixLB_internal <- function(
	### internal function called on each cluster node from \code{\link{twSfApply2MatrixLB}}
	index,XNamesNames,FUN,...
){
	X1 <- eval.parent(XNamesNames[[1]])
	X2 <- eval.parent(XNamesNames[[2]])
	FUN(X1[index,],X2[index,],...)
}


