#---------- process knowledge: the model relating microbial parameters to observations

modMetaICBM1 <- function(
### Creating Meta-information for ICBM1.
){
	##seealso<< \code{\link{twCreateModMeta}}, \code{\link{twSetModMetaAuxGroups}}
	modMeta <- twCreateModMeta(
		rowNames =  c('Y','O','R')
		#,csts = list( cis=c('c12','c14'), nis = c('n15') )
		,csts = list( cis=c('c12','c14') )
		,iRUnits=c(c12=1,c14=c14Constants$iR14CUnit)		
	)
	modMeta$iR14CStandard <- c14Constants$iR14CStandard
	
	auxGroups <- list(
		#csums = paste("csums",modMeta$rowNames,sep="_")	#sum across all c compartments 
		inputLeaf = paste("inputLeaf",modMeta$colNames,sep="_")
		,inputRoot = paste("inputRoot",modMeta$colNames,sep="_")
		,decY = paste("decY",modMeta$colNames,sep="_") 
		,decO = paste("decO",modMeta$colNames,sep="_")
		,respY = paste("respY",modMeta$csts$c,sep="_") 			
		,respO = paste("respO",modMeta$csts$c,sep="_") 	
		#,dxSums = paste("dxSums",modMeta$csts$c,sep="_")	#sum across all c compartments
	)
	auxGroupsSolve <- list(
		cStock = "cStock"								#sum over 12C of all SOM pools 
		,F14C = paste("F14C",modMeta$rowNames,sep="_")	#fraction modern of soil carbon pools 			
		,F14CT = "F14CT"					#fraction modern of entire soil (mixing model across pools) 			
		#,respF14C = paste("respF14C",modMeta$rowNames,sep="_") 	#fraction modern of respiration from soil carbon pools
		,respF14CT = "respF14CT" 					#fraction modern of total respiration
	)
	modMeta <- twSetModMetaAuxGroups(modMeta,auxGroups)
	#copy2clip(paste("enum AUX_OUTPUT_NAMES {",paste(c(modMetaICBM1()$auxOutputNames,"N_AUX"),collapse=","),"}; //generated in ICBM1.R",sep=""))	#to adjust icbm1.c 
}
#twUtestF("ICBM1",test="modMeta")
#mtrace(modMetaICBM1)

initStateICBM1 <- function(
	### Creating initial state variables for ICBM1.
	xc12,	cn=0, 	iR 
	,modMeta=modMetaICBM1()	##<< may pass pre-calulated modMeta for efficiency.
	,... 
){
	##seealso<< 
	x <- initStateSoilMod(xc12,cn,iR,modMeta=modMeta)
	### Numeric matrix (nPool, nIsotopes) of state variable mass.
}
#twUtestF("ICBM1",test="init")

#mtrace(solveICBM1)
solveICBM1 <- function(
	### solve the ODE of \code{\link{derivICBM1}}
	x0		##<< numeric vector or matrix at t=0
	,times	##<< times at which explicit estimates for y are desired. The first value in times must be the initial time.
	,parms	##<< list of model parameters
	,input  ##<< list with dataframes entries leaf and root each with columns yr and obs
	,fFmAtmosphere=fmAtmosphere  
	,modMeta=modMetaICBM1()	##<< metaInformation from model. Pass for efficiency or when using different units. 
	,useRImpl=FALSE		##<< flag indicating to use the R implementation instead of C implementation.
){
	##seealso<< \code{\link[deSolve]{lsoda}}
	# Tried to works with lsoda(..., atol=0). Else microbial biomass may become truely zero insteald of
	# small values, and derivative functions produces NaNs. But still gets 0
	parms$modMeta <- modMeta	#a way to pass it to the derivative function
	# -- create the forcing according to time lag for root input 14C
	fLeaf12C <- if( nrow(input$leaf) == 1) function(t){as.vector(input$leaf[1,2])} else
	   approxfun(x=input$leaf[,1], y=input$leaf[,2], method="linear", rule=2)
	fRoot12C <- if( nrow(input$root) == 1) function(t){as.vector(input$root[1,2])} else
	   approxfun(x=input$root[,1], y=input$root[,2], method="linear", rule=2)
	#for 14C we need an input each year, because it changes with atmospheric 14c-Activity
	iTimes <- seqRange( range(times), by=1 )
	input$leaf14C <- cbind( yr=iTimes, obs=fLeaf12C(iTimes)*fFmAtmosphere(iTimes-parms$tLagLeaf) ) 
	input$root14C <- cbind( yr=iTimes, obs=fRoot12C(iTimes)*fFmAtmosphere(iTimes-parms$tLagRoot) ) 
	res0 <- if( useRImpl ){ 
		#lsoda( x0, times, derivICBM1, parms, atol = 0 )
		# the forcing functions; rule = 2 avoids NaNs in interpolation
		fLeaf14C = approxfun(x=input$leaf14C[,1], y=input$leaf14C[,2], method="linear", rule=2)
		fRoot14C = approxfun(x=input$root14C[,1], y=input$root14C[,2], method="linear", rule=2)
		parms$fInput <- function(t){ 
					c(leaf12C = fLeaf12C(t), leaf14C=fLeaf14C(t), root12C=fRoot12C(t), root14C=fRoot14C(t))    
				}
		lsoda( x0, times, derivICBM1, parms ) 
	}else {
		forcings <- input[c("leaf","root","leaf14C","root14C")] 
		lsoda( x0, times, dllname = "howlandInversion", func = "deriv_icbm1",	initfunc = "init_soilmod_icbm1", parms = parms, nout = modMeta$nAux, outnames = modMeta$auxOutputNames, initforc = "forcc_icbm1", forcings=forcings)		#head(res0 <- lsoda( x0, times, dllname = "howlandInversion", func = "deriv_icbm1",	initfunc = "init_soilmod_icbm1", parms = parms, nout = modMeta$nAux, outnames = modMeta$auxOutputNames, initforc = "forcc_icbm1", forcings=forcings))
	}
	#cname=modMeta$rowNames[1]
	auxF14C <- do.call( cbind, lapply( modMeta$rowNames, function(cname){(res0[,paste(cname,"c14",sep="_")]/res0[,paste(cname,"c12",sep="_")]) } )) / modMeta$iR14CStandard
	colnames(auxF14C) <- paste("F14C",modMeta$rowNames,sep="_")
	cStock <- rowSums(res0[,c("Y_c12","O_c12")])
	res <- cbind( res0 
			#,inputLeaf12C=leaf12C(times), inputLeaf14C=leaf14C(times), inputRoot12C=root12C(times), inputRoot14C=root14C(times)
			,cStock=cStock
			,F14CT=rowSums(auxF14C[,c("F14C_Y","F14C_O")]*res0[,paste(c("Y","O"),c("c12"),sep="_")]) / cStock 
			,respF14CT=rowSums(res0[,c("respY_c14","respO_c14")])/rowSums(res0[,c("respY_c12","respO_c12")]) / modMeta$iR14CStandard
			,auxF14C
		)
	### result of \code{\link{lsoda}}  
}

derivICBM1 <- function(
	### Derivative function of Basic Colimitation model.
	t, x, p 
){
	##details<< 
	## Calls forcing function p$fInput which should (see \code{\link{forcings}}). 
	# Meta information of this model passed with p
	mm <- p$modMeta
	# the derivative
	dx <- mm$matrixTemplate
	# auxiliary outputs 
	a <- mm$auxOutputTemplate
	
	# ode solvers provide x as a vector, attach attributes to have access as a mtrix.
	if( !is.matrix(x) ){ 
		dim(x) <- c(mm$nRow, mm$nCol)
		dimnames(x) <- dimnames(dx)
	}
	
	decY <- p$kY * x["Y",]
	respY <- decY[mm$csts$cis] * (1-p$h)  
	decO <- p$kO * x["O",]
	respO <- decO[mm$csts$cis]
	resp <- rbind(respY,respO)
	input <- p$fInput(t)
	inputLeaf <-  c(input["leaf12C"],input["leaf14C"])  
	inputRoot <-  c(input["root12C"],input["root14C"])  
	dx["Y",] <-  +inputLeaf +inputRoot - decY  
	dx["O",] <-  +p$h*decY - decO
	dx[,"c14"] <- dx[,"c14"] -x[,"c14"]*c14Constants$lambda	#radioactive decay see c14Utils
	dx["R",] <- +respY +respO
	
	#a[ mm$auxGroups$csums ] = csums
	a[ mm$auxGroups$inputLeaf ] = inputLeaf
	a[ mm$auxGroups$inputRoot ] = inputRoot
	a[ mm$auxGroups$decY ] = decY
	a[ mm$auxGroups$decO ] = decO
	a[ mm$auxGroups$respY ] = respY
	a[ mm$auxGroups$respO ] = respO
	#a[ mm$auxGroups$cStock] = sum(x[c("Y","O"),"c12"])
	#a[ mm$auxGroups$F14C] = (x[,"c14"]/x[,"c12"]) / mm$iR14CStandard
	#a[ mm$auxGroups$F14CT] = sum(a[mm$auxGroups$F14C][c("F14C_Y","F14C_O")]*x[c("Y","O"),"c12"])/a[ mm$auxGroups$cStock]  
	#a[ mm$auxGroups$respF14C] = (resp[,"c14"]/resp[,"c12"]) / mm$iR14CStandard
	#a[ mm$auxGroups$respF14CT] = (sum(resp[,"c14"])/sum(resp[,"c12"])) / mm$iR14CStandard
	
	list(dx,a)
}

calcSteadyK_ICBM1 <- function(
	### calculte decay constants from assuming steady state and remaining parameters
	Ctot	##<< SOM C-Stock
	,iY		##<< steady state input
	,parms	##<< list with entries kY,kO and cY and h
){
	padj <- parms
	padj$kY = iY/(padj$cY*Ctot)
	padj$kO = (iY*padj$h)/((1-padj$cY)*Ctot)
	padj
	### named numeric matrix with columns kY and kO, rows corresponding to longest input parameter (others are recycled)
}

calcSteadyK_ICBM1_mat <- function(
	### calculte decay constants from assuming steady state and remaining parameters
	Ctot	##<< SOM C-Stock
	,iY		##<< steady state input
	,parms	##<< named numeric vector with entries kY,kO and cY and h
){
	padj <- parms
	padj[,"kY"] = iY/(padj[,"cY"]*Ctot)
	padj[,"kO"] = (iY*padj[,"h"])/((1-padj[,"cY"])*Ctot)
	padj
	### named numeric matrix with columns kY and kO, rows corresponding to longest input parameter (others are recycled)
}


calcSteadyHcY_ICBM1 <- function(
	### calculte cY and h from assuming steady state and decay constants
	Ctot	##<< SOM C-Stock
	,iY		##<< steady state input
	,parms	##<< named numeric vector with entries kY,kO and cY and h
){
	padj <- parms
	padj$cY <- cY <- iY/padj$kY/Ctot
	padj$h <- padj$kO*(1-cY) / (padj$kY*cY)
	padj$Ctot0 <- Ctot
	padj
	### parms with entries cY and h replaced
}

calcSteadyHcY_ICBM1_mat <- function(
	### calculte cY and h from assuming steady state and decay constants, set Ctot0 to Ctot
	Ctot	##<< SOM C-Stock
	,iY		##<< steady state input
	,parms	##<< named numeric matrix with columns kY,kO and cY and h
){
	padj <- parms
	padj[,"cY"] <- cY = iY/padj[,"kY"]/Ctot
	padj[,"h"] <- padj[,"kO"]*(1-cY) / (padj[,"kY"]*cY)
	padj[,"Ctot0"] <- Ctot 
	padj
	### parms with columns cY,h,Ctot0 replaced
}

calcRelaxSteadyHcY_ICBM1 <- function(
	### calculate cY,h, from assuming increase rate dO of the old pool and decay constants and given C0
	Ctot	##<< not used (uses parms$Ctot0)
	,iY		##<< sum of inputs
	,parms	##<< named numeric vector with entries kY,kO,dO,Ctot and cY and h
){
	##details<<
	## The young pool is assumed to be in steady state.
	## The old pool is assumed to increase with rate h
	## When calculating litter input, then roughly input = dO + respiration
	## This parameter calculation function does calculate Ctot0. Make sure that Ctot0 is among estimated parameters. Else initialization function will fail. 
	padj <- parms
	Ctot = parms$Ctot0
	padj$cY = cY = iY/padj$kY/Ctot
	padj$h = (padj$dO/Ctot + padj$kO*(1-cY)) / (padj$kY*cY)
	padj
	### parms with entries cY and h replaced
}

calcRelaxSteadyHcYC0_ICBM1 <- function(
	### calculate cY,h, and initial Ctot from assuming increase rate dO of the old pool and decay constants and approximate C0
	Ctot	##<< SOM C-Stock in yr parms$yrCtot
	,iY		##<< here the respiration, dO will be added
	,parms	##<< named numeric vector with entries kY,kO,dO,yr0 (initial time),yrCtot (time of Ctot measurement) and cY and h
){
	##details<<
	## Ctot0 is calculated as Ctot - (yrCtot-yr0)*dO
	## the results from \code{\link{calcRelaxSteadyHcY_ICBM1}} are returned
	padj <- parms
	dt <- padj$yrCtot - padj$yr0
	padj$Ctot0 <- Ctot0 <- max(1e-4, Ctot - padj$dO*dt)
	calcRelaxSteadyHcY_ICBM1( NA, iY=iY, parms=padj )
	### parms with entries cY,h, and Ctot replaced
}

calcSteadyCtot0CY_ICBM1 <- function(
	### calculate Ctot0 and cY from assuming yong pool in steady state and C0 from linear change to Ctot
	Ctot	##<< SOM C-Stock
	,iY		##<< steady state input
	,parms	##<< named numeric vector with entries kY,kO,dO,yr0,yrCtot 
){
	##details<< 
	## cY is the proportion Y/C0 at the time of yr0
	##
	## 
	padj <- parms
	deltaT <- padj$yrCtot - padj$yr0
	padj$Ctot0 <- Ctot0  <- max(1e-8, Ctot - deltaT*padj$dO)
	padj$cY <- cY <- iY/padj$kY/Ctot0
	padj
	### parms with entries C0 and cY replaced
}










