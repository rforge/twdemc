data(dsKinrespTut)
#library(nlme)
#source(file.path("R","kinrespModel1.R"))

#---------  the data: respiration measurements at soil core scale
dsKinrespTut[1:10,]		# first 10 records
summary(dsKinrespTut)	# variables and overview their distributions

windows(width=3.4,height=2.8, pointsize=10, record=TRUE)
	par( las=1 )					#also y axis labels horizontal
	par(mar=c(2.2,2.3,0,0)+0.3 )  #margins
	par(tck=0.02 )				#axe-tick length inside plots             
	par(mgp=c(1.4,0.2,0) )  #positioning of axis title, axis labels, axis
plot( resp ~ time, data=dsKinrespTut, col=dsKinrespTut$replicate, xlab="Time (hours)", ylab="Respiration (µg Resp-C/g Soil/hour)", pch=as.numeric(experiment) )


ds <- subset(dsKinrespTut, experiment==9 & replicate==3 )

#--------- the model
?kinrespModel

#------------ least squares fit to one replicate from package,
# for fitting several replicates use nlme including random effects
# for simplicity we assume iid errors here

#?gnls		# help on the gnls function
fitLS <- gnls( resp ~ kinrespModel(x0,r0,mumax,time), ds
	,params=x0+r0+mumax~1
	,start=c(x0=140, r0=0.1, mumax=0.24) )
coef(fitLS)


#fitLS2$par			# compare the LS-Hand estimates
confint(fitLS)		# from LogLik profile assming gaussian error distribution
summary(fitLS)

# plot the data and the fit together
plot( resp ~ time, data=ds, col=replicate, xlab="Time (hours)", ylab="Respiration (µg Resp-C/g Soil/hour)", pch=as.numeric(experiment) )
lines( fitted(fitLS)~time, data=ds )

# check model assumptions of iid residuals
# there is a problem
plot((resid(fitLS))~ds$time)	# high leverage of last observations, magnitude of error increases
plot((resid(fitLS)/fitted(fitLS))~ds$time)	# relative error decreases

#-------------- Additive generalized model: error, increase with magnitude of flux
fitLS5 <- gnls( resp ~ kinrespModel(x0,r0,mumax,time), ds
	,params=x0+r0+mumax~1
	,start=c(x0=140, r0=0.1, mumax=0.24)
	,weights=varPower(fixed=0.5)
)
coef(fitLS5)
coef(fitLS)			# compare (not very different)
plot((resid(fitLS5,type="pearson"))~ds$time)	# seems ok for weighted residuals

# plot the data and the fit together
plot( resp ~ time, data=ds, col=replicate, xlab="Time (hours)", ylab="Respiration (µg Resp-C/g Soil/hour)", pch=as.numeric(experiment) )
lines( fitted(fitLS5)~time, data=ds, col="maroon" )
