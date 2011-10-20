findSplit <- function(
	### determine the parameter and its value, where to split the parameter space
	ss	##<< the sample (stacked chains)
	,nSplit = 4 ##<< how points per parameter dimension to check for split
	,rVarCrit = 20^2	##<< minimla ratio of variances of a variable between right and left side of a splitting point
	,rAlphaSlopeCrit = base:::pi/2/3	## minimal ratio of angle between scaled slopes, defaults to a third of a quarter of a cirle 
	,iVars = 1:ncol(ss)		#<< variables to check for split		
	,jVars = 1:ncol(ss)		#<< variables to check for different scales of variance
){
	##details<< 
	## First it checks for different scales of variance in other variables. 
	foundSplit <- FALSE
	quantVar <- iVarLeft <- iVarRight <- matrix( NA_real_, nrow=length(iVars), ncol=nSplit )
	slope <- parVarLeft <- parVarRight <- array( NA_real_, dim=c(length(iVars),length(jVars),nSplit)
	, dimnames=list(iVar=colnames(ss)[iVars], jVar=colnames(ss)[jVars], iSub=NULL ))
	ssSubLeft <- ssSubRight <- vector(mode="list",length=length(iVars) )
	#i <- 1
	for( i in seq_along(iVars)){
		iVar <- iVars[i]
		p <- ss[,iVar] 
		qp <- quantVar[i, ] <- quantile( p, 1:(nSplit)/(nSplit+1) )
		#iSub <- 1
		#iSub <- nSub
		ssSubLeft[[i]] <- lapply( 1:(nSplit), function(iSplit){
				ssSub <- ss[ (p <= qp[iSplit]),]
			})
		ssSubRight[[i]] <- lapply( 1:(nSplit), function(iSplit){
				ssSub <- ss[ (p >= qp[iSplit]),]
			})
		# ratio of variances of subsets left and right of splitting points
		#j <- 1
		for( j in seq_along(jVars)){
			jVar <- jVars[j]
			#iSplit=1
			if( jVar!=iVar){
				pVar <- rep( NA_real_, nSplit)
				for(iSplit in 1:nSplit){
					varLeft <- parVarLeft[i,j,iSplit] <- var(ssSubLeft[[i]][[iSplit]][,jVar])
					varRight <- parVarRight[i,j,iSplit] <- var(ssSubRight[[i]][[iSplit]][,jVar]) 
					#plot( ssSubLeft[[i]][[iSplit]][,iVar], ssSubLeft[[i]][[iSplit]][,jVar] )
					pVar[iSplit] <- varLeft /varRight
				}
				pVar[ pVar < 1 ] <- 1/pVar[pVar<1]
				if( max(pVar) > rVarCrit){
					iSplit <- which.max(pVar)
					return( structure( quantVar[i,iSplit], names=colnames(ss)[iVar] ) )
				}
			} # end if( jVar!=iVar)
		} # end for jVar variance
	} # for iVar
	
	##details<< 
	## Next it checks for different angles of the normalized slopes in relation of the parameters 
	#i<-1
	if( !foundSplit) for( i in seq_along(iVars)){
		# check difference in slope angles
		iVar <- iVars[i]
		# also need the variance of iVar for scaling the slope
		iPos <- match( iVar, jVars )
		iVarLeft[i, ] <- sapply( 1:(nSplit), function(iSplit){ var(ssSubLeft[[i]][[iSplit]][,iVar]) })
		iVarRight[i, ] <- sapply( 1:(nSplit), function(iSplit){ var(ssSubRight[[i]][[iSplit]][,iVar]) })
		#qp <- quantVar[i, ]
		#j <- 1
		for( j in seq_along(jVars)){
			jVar <- jVars[j]
			if( jVar!=iVar){
				#iSplit=4
				dAlphaSlope <- sapply( 1:nSplit, function(iSplit){
						#XXTODO names are not know to construct formula, lookup direct equation for
						#plot( ssSubLeft[[i]][[iSplit]][,iVar], ssSubLeft[[i]][[iSplit]][,jVar] )
						lmLeft <- lm( as.formula(paste(colnames(ss)[c(jVar,iVar)],collapse=" ~ ")), data=as.data.frame(ssSubLeft[[i]][[iSplit]][,c(jVar,iVar)]) )
						lmRight <- lm( as.formula(paste(colnames(ss)[c(jVar,iVar)],collapse=" ~ ")), data=as.data.frame(ssSubRight[[i]][[iSplit]][,c(jVar,iVar)]) )
						#plot( ssSubLeft[[i]][[iSplit]][,iVar], ssSubLeft[[i]][[iSplit]][,jVar] ); abline(lmLeft)
						#plot( ssSubRight[[i]][[iSplit]][,iVar], ssSubRight[[i]][[iSplit]][,jVar] ); abline(lmRight)
						slopeLeft <- coef(lmLeft)[2] * sqrt(iVarLeft[i,iSplit]/parVarLeft[i,j,iSplit]) 
						slopeRight <- coef(lmRight)[2] * sqrt(iVarRight[i,iSplit]/parVarRight[i,j,iSplit])
						alphaRight <- asin( slopeLeft)
						alphaLeft <- asin( slopeRight )
						abs( alphaLeft - alphaRight )
					})
				if( max(dAlphaSlope) > rAlphaSlopeCrit){
					iSplit <- which.max(dAlphaSlope)
					return( structure( quantVar[i,iSplit], names=colnames(ss)[iVar] ) )
				}
			} # end if( jVar!=iVar)
		} # end for jVar variance
	} # for iVar

	### named scalar: splitting value with the name of the parameter dimension that is to split 
	return( NA_real_ )
}
attr(findSplit,"ex") <- function(){
	ss1 <- ss <- pps[,-1]
	#mtrace(findSplit)
	(res <- findSplit(ss1))
		
	ss2 <- ss <- ss1[ ss1[,res$iVar] > res$splitValue,]
	res2 <- findSplit(ss2)
	ss3 <- ss <- ss2[ ss2[,res2$iVar] > res2$splitValue,]
	res3 <- findSplit(ss3)
	ss4 <- ss <- ss3[ ss3[,res3$iVar] > res3$splitValue,]
	res4 <- findSplit(ss4)
	
	#mtrace(findSplit)
	(res <- findSplit(ss1, iVars=c(2,1) ))
	(res <- findSplit(ss1, iVars=c(2) ))
	plot( ss[,1], ss[,2], col=c("blue","red")[as.integer(ss[,names(res)] >= res)+1])
}