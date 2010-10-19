# AP = analytical precision
typical_precision <- 0.005	# analytical
high_precision <- 0.001		# analytical
precisionSpatHet <- 0.01	# including spatial heterogeneity
precision <- precisionSpatHet


#==============================================================================
# Howland data 
#==============================================================================
# Howland site, values measured with ruler from 2-pool graph from Howland2009summary_fromSueTrumbore.doc

# unit in gC/m2 and yr
### Observations at Howland 
Howland14C <- list(
	obsNutrientSite = list( 
		##describe<<
		respCum = cbind(times = 2007, obs = 365, sdObs = 30)	##<< cumulated carbon in respiration from Reservour-C (non-recent C) in gC/m2
		,respFM = {tmp=delta2FM(127.5); cbind(times=1999, obs = tmp, sdObs = precision*tmp)} ##<< fraction modern of respiration	
		,somStock = cbind(times = 2007, obs = 1010, sdObs = 100) #AP not given!	##<< C-stock in SOM gC/m2
		,somOStock = cbind(times = 2007, obs = 920, sdObs = 100) #AP not given!	##<< C-stock in O-horizon gC/m2
		,somOFM = {tmp=delta2FM(c(225,195,150)); cbind(times = c(1988,1996,2007), obs = tmp , sdObs = precision*tmp)} ##<< fraction modern of SOM (O-Horizon)
		##end
	)
	,obsTowerSite = list( 
		##describe<<
		respCum = cbind(times = 2007, obs = 365, sdObs = 30)	##<< cumulated carbon in respiration from Reservour-C (non-recent C) in gC/m2
		,respFM = {tmp=delta2FM(150); cbind(times=1999, obs = tmp, sdObs = precision*tmp)} ##<< fraction modern of respiration	
		,somStock = cbind(times = 2007, obs = 10100, sdObs = 100) #AP not given!	##<< C-stock in SOM gC/m2
		,somOStock = cbind(times = 2007, obs = 6500, sdObs = 100) #AP not given!	##<< C-stock in O-horizon gC/m2
		,somOFM = {tmp=delta2FM(c(120,118)); cbind(times = c(1996,2007), obs = tmp, sdObs = precision*tmp)} ##<< fraction modern of SOM (O-Horizon)
		##end
	)
	,litter = list(
		leaf = cbind(times=1997, obs=230, sdObs=40)
		,root = cbind(times=1997, obs=250, sdObs=100)
	)
)
save(Howland14C,file=file.path("data","Howland14C.RData"))


