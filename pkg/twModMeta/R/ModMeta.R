twCreateModMeta <- function(
	### Creating meta-information about rows (stateVars), isotopes (columns), and additonal outputs. 
	rowNames, 		##<< character vector: names of the state variables 
	csts, 	  		##<< list of character vectors: correponding to isotopes
	iRUnits=1, 	  	##<< numeric vector of units for different isotopes, ensuring equal magnitudes to avoid numerical errors (see details)\cr
	auxGroups=list() ##<< namegroup (character) -> columnNames (character vector): column names of auxiliary outputs 
){
	if( !is.vector(rowNames) ) 	stop("rowNames must be a vector of character strings of at least length one.")
	if( !(is.character(rowNames)) ) stop("rowNames must be a vector of character strings of at least length one.")
	if( !(is.list(csts)) ) stop("csts must be a list of vectors of character strings of at least length one.")
	#if( !(is.character(csts[[1]][1])) ) stop("csts must be a list of vectors of character strings of at least length one.")
	res <- within( list(),{
			rowNames <- as.character(rowNames)	#strip names attribute	
			csts <- csts
			colNames <- as.character(unlist(csts))
			nRow <- length(rowNames)
			nCol <- length(colNames)
			elementNames <- twModElementNames(rowNames,colNames)
			matrixTemplate <- matrix( 0.0, nrow=nRow, ncol=nCol, dimnames=list(pools=rowNames, compartments=colNames) )
		})
	names(res$matrixTemplate) <- res$elementNames
	##details<< \describe{\item{iRUnits}{ 
	## By default all units are 1.
	## If names are not given, is is assumed that order corresponds to colNames
	## If names are given, default 1 for c12 and n15 are added.
	## }}
	if( 0<length(iRUnits) ) iRUnits=1
	if( length(iRUnits)==1 ) 
		iRUnits <- structure( rep(iRUnits,res$nCol), names=res$colNames )
	else{
		if( is.null(names(iRUnits))){ 
			if( length(iRUnits != res$nCol) ) stop("iRUnits must be provided for each isotope component")
			names(iRUnits) <- res$colNames	#assume that order is correct
		}else{	
			if( !("c12" %in% names(iRUnits)) ) iRUnits <- c(iRUnits,c12=1)	#add default entry for c12
			if( !("n15" %in% names(iRUnits)) ) iRUnits <- c(iRUnits,n15=1)	#add default entry for n15
			iRUnits <- iRUnits[res$colNames]	#permute to correct order 
			if( length(iRUnits != res$nCol) ) stop("iRUnits must be provided for each isotope component")
		}
	}
	res$iRUnits <- iRUnits	
	res <- twSetModMetaAuxGroups(res, auxGroups)
	### a list with components \describe{
	### \item{rowNames}{character vector: names of the state variables }
	### \item{csts}{ list of character vectors: correponding to isotopes }
	### \item{colNames}{character vector: names of the isotopes}
	### \item{elementNames}{character vector: names of the elements row_column }
	### \item{iRUnits}{numeric array (nCol), unit of isotopic column, take care in sums}	
	### \item{auxGroups}{namegroup (character) -> columnNames (character vector): column names of auxiliary outputs }
	### \item{auxOutputNames}{character vector: names of all auxiliary outputs}	
	### \item{auxOutputTemplate}{numeric vector, with names corresponding to auxOutputNames }	
	### \item{nRow}{number of rows}	
	### \item{nCol}{number of columns}	
	### \item{matrixTemplate}{numeric matrix: template for state variable matrix with corresponding row and column names}	
	### \item{nAux}{number of auxiliary outputs}	
	### }
}

twSetModMetaAuxGroups <- function(
	### Setting auxiliare output items in modMeta.
	modMeta, auxGroupsNew){
	res <- within(modMeta,{
			auxGroups <- auxGroupsNew			#grouping names of auxiliary outputs for assigning vectors
			auxOutputNames <- as.character(unlist(auxGroups))
			if( 0 == length(auxOutputNames) ) auxOutputNames <- character()
			auxOutputTemplate <- structure( rep(0.0,length(auxOutputNames)), names=auxOutputNames)
			nAux <- length(auxOutputNames)
		})
	### ModMeta (result of \code{\link{twCreateModMeta}}) adjusted for auxGroups, ausOutputnames, auxOutputTemplate, and nAux.
}
#mtrace(twSetModMetaAuxGroups)

twModElementNames <- function(
	### Construct element names for given rows and columns
	rowNames, ##<< A character vector of row names.
	colNames  ##<< A character vector of column names.
){ 
	as.vector(sapply( colNames, function(colName){ paste(rowNames, colName, sep="_" )}))
	### Character vector for each index (by row) rowName, colName
}

DLLfuncTest <- function (
	### Adjustment of \code{\link[deSolve]{DLLfunc}}, that does not coerce parms to double vector.
	func, times, y, parms, dllname, initfunc = dllname, 
	rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL, 
	initforc = NULL, fcontrol = NULL) 
{
	if (!is.numeric(y)) 
		stop("`y' must be numeric")
	n <- length(y)
	if (!is.null(times) && !is.numeric(times)) 
		stop("`times' must be NULL or numeric")
	if (!is.null(outnames)) 
		if (length(outnames) != nout) 
			stop("length outnames should be = nout")
	ModelInit <- NULL
	Outinit <- NULL
	flist <- list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
	Ynames <- attr(y, "names")
	if (is.null(dllname) || !is.character(dllname)) 
		stop("`dllname' must be a name referring to a dll")
	if (!is.null(initfunc)) 
		if (is.loaded(initfunc, PACKAGE = dllname, type = "") || 
			is.loaded(initfunc, PACKAGE = dllname, type = "Fortran")) {
			ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
		}
		else if (initfunc != dllname && !is.null(initfunc)) 
			stop(paste("cannot integrate: initfunc not loaded ", 
					initfunc))
	if (is.null(initfunc)) 
		initfunc <- NA
	#if (!is.null(forcings)) 
#		flist <- checkforcings(forcings, times, dllname, initforc, 
#			TRUE, fcontrol)
	if (!is.character(func)) 
		stop("`func' must be a *name* referring to a function in a dll")
	if (is.loaded(func, PACKAGE = dllname)) {
		Func <- getNativeSymbolInfo(func, PACKAGE = dllname)$address
	}
	else stop(paste("cannot run DLLfunc: dyn function not loaded: ", 
				func))
	dy <- rep(0, n)
	storage.mode(y) <- storage.mode(dy) <- "double"
	out <- .Call("call_DLL", y, dy, as.double(times[1]), Func, 
		ModelInit, parms, as.integer(nout), as.double(rpar), 
		as.integer(ipar), as.integer(1), flist, PACKAGE = "deSolve")
	vout <- if (nout > 0) 
			out[(n + 1):(n + nout)]
		else NA
	out <- list(dy = out[1:n], var = vout)
	if (!is.null(Ynames)) 
		names(out$dy) <- Ynames
	if (!is.null(outnames)) 
		names(out$var) <- outnames
	return(out)
}

initStateModMeta <- function(
	### Creating initial state variables for Basic Colimitation submodel.
	xc12 		
	### numeric vector: 12C mass for each pool \cr
	### If names are given these must comprise code{modMeta$rowNames}.
	,cn 		
	### numeric vector: carbon to nitrogen ratio for each pool \cr
	### If vector has length one, the same ratio is assumed for all pools.
	### If names are given, these must comprise code{modMeta$rowNames}
	,iR 
	### numeric matrix (nPools, nIsotopes): atomic ratios \cr
	### For others be aware of the adjusted units to avoid bigger numerical errors.
	### If it is a vector (nIsotopes), the same atomic ratio is assumed for all pools   
	### If colnames/names are given they comprise code{modMeta$colNames}
	### Default 1 is added for columns c12 and n15.
	,modMeta	##<< model meta information, see \code{\link{modMetaICBMDemo}}
){
	# check arguments
	#rowNamesSOM <- modMeta$rowNames[ modMeta$rowNames != "R"]
	rowNamesSOM <- modMeta$rowNames
	if( !is.null(names(xc12))) xc12 <- xc12[rowNamesSOM]	#permute items of xc12	
	if( (modMeta$nRow) != length(xc12) ) stop("xc12 must contain a value for each state variable.")
	if( !is.matrix(iR) )	#repeat for each row
		iR <- matrix( iR, byrow=TRUE, nrow=modMeta$nRow, ncol=length(iR), dimnames=list(rowNamesSOM,names(iR)))
	if( !is.null(rownames(iR))) iR <- iR[rowNamesSOM,,drop=FALSE]	#permute rows 
	if( (modMeta$nRow) != nrow(iR) ) stop("iR must contain a row for each state variable.")
	
	if( !("c12" %in% colnames(iR)) ) iR <- cbind(iR,c12=1)	#add default entry for c12
	iRcis <- iR[,modMeta$csts$cis,drop=FALSE]				#extract carbon and permute 
	if( length(modMeta$csts$cis) != ncol(iRcis) ) stop("iR must contain a value for each carbon isotope.")
	if( !("n15" %in% colnames(iR)) ) iR <- cbind(iR,n15=1)
	iRnis <- iR[,modMeta$csts$nis,drop=FALSE]
	if( length(modMeta$csts$nis) != ncol(iRnis) ) stop("iR must contain a value for each nitrogen isotope.")
	
	if( length(cn) == 1)		#repeat for each row
		cn <- rep(cn,length(rowNamesSOM))
	if( !is.null(names(cn)) ) cn <- cn[rowNamesSOM]	#permute names
	if( !length(cn) == (modMeta$nRow)) stop("cn must contain a value for each state variable.")
	
	#initialize the state variables
	x <- modMeta$matrixTemplate
	x[ rowNamesSOM ,modMeta$csts$cis] <- xc12*iRcis
	xn15 <- xc12*cn
	x[ rowNamesSOM, modMeta$csts$nis] <- xn15*iRnis
	x
	### Numeric matrix (nPool, nIsotopes) of state variable mass.
}

twStateMatODERes <- function(
	### returns the matrix of state variables for an ode output corresponding to modMeta
	out			##<< result of ode to model of matrix statevariables
	,iSteps		##<< integer vector index of outputs
	,modMeta	##<< description of the model
){
	#iSteps=1
	#iSteps=1:3
	#iStep=1
	res <- abind( lapply( iSteps, function(iStep){
		matrix( out[iStep,2:(modMeta$nCol*modMeta$nRow+1)], nrow=modMeta$nRow )
	}), rev.along=0 )
	dimnames(res) <- c( modMeta[c("rowNames","colNames")], list(outStep=NULL) )
	res
	### array (nRow,nCol,length(iSteps)) of stacked state variable matrices for each step 
}

