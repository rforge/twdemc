data(delta14Catm)
#source(file.path("R","c14Utils"))

### Constants used with 14C calculations.
c14Constants <- list(
	#iR14CUnit = 1e-12
	#,iR14CStandard = 1.176	# http://www.nosams.whoi.edu/clients/data.html,  (Karlen, et. al., 1968),
	iR14CUnit = 1.176*1e-12		##<< atomic ratio of standard, i.e. output in fraction modern
	,iR14CStandard = 1			# http://www.nosams.whoi.edu/clients/data.html,  (Karlen, et. al., 1968), 
	,yr14CStandard=1950
	,lambda=1/8267	#radioactive decay constant
	,delta14Catm=delta14Catm
	,fmAtm=cbind(yr=delta14Catm[,"yr"], fm14C=delta2FM(delta14Catm[,"delta14C"]))	##<< fraction modern of the athmosphere
)

save(c14Constants, fmAtmosphere, file=file.path("data","c14Constants.RData"))

