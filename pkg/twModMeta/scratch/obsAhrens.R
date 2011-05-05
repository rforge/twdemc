# AP = analytical precision
typical_precision <- 0.005
high_precision <- 0.001
precision <- high_precision

#==============================================================================
# Howland data 
#==============================================================================
# Howland site, values measured with ruler from 2-pool graph from Howland2009summary_fromSueTrumbore.doc

#------------------------------------------------------------------------------
# Nutrient control site
#------------------------------------------------------------------------------
obsFM_CO2 <- cbind(time = 1999, FM_CO2 = 127.5/1000 + 1, AP = precision)
obsFM_SOM <- cbind(time = c(1988,1996,2007), FM_SOM = c(225,195,150)/1000 + 1, AP = precision)
obsC_tot <- cbind(time = 2007, C_tot = 1010, AP = 100) #AP not given!
obsRes_CO2 <- cbind(time = 2007, Res_CO2 = 365, AP = 30) 
# from Howland2009_summary. page 2

#------------------------------------------------------------------------------
# Tower control site
#------------------------------------------------------------------------------
#obsFM_CO2 <- cbind(time = 1999, FM_CO2 = 150/1000 + 1, AP = precision)
#obsFM_SOM <- cbind(time = c(1996,2007), FM_SOM = c(120,118())/1000 +1, AP = precision)
#obsC_tot <- cbind(time = 2007, C_tot = 10100, AP = 100) #AP not given!
#
## from Howland2009_summary. page 2
#obsRes_CO2 <- cbind(time = 2007, Res_CO2 = 365, AP = 30) 


## from Davidson et al. 2006: Table 1, ResCO2 [gC m-2] = 100*0.5*soil respiration [MgC ha-1]
#obsRes_CO2 <- cbind(time = c(1997,1998,1999,2000,2001,2002,2003), Res_CO2 = c(6.1,7.6,7.5,7.3,7.3,7.5,8.1)*100*0.5, AP = 30)  


# combining all observation in a complete time series where gaps are NAs
times <- seq(1900, 2007, by = 1)

obs <- cbind(time = times, obsFM_CO2 = NA, obsFM_SOM = NA, obsC_tot = NA, obsRes_CO2 = NA) 
selection <- obs[,"time"] %in% obsFM_CO2[,"time"]

obs[selection,2] <- obsFM_CO2[,"FM_CO2"]
obs[obs[,"time"] %in% obsFM_SOM[,"time"],"obsFM_SOM"] <- obsFM_SOM[,"FM_SOM"]
obs[obs[,"time"] %in% obsC_tot[,"time"],"obsC_tot"] <- obsC_tot[,"C_tot"]
obs[obs[,"time"] %in% obsRes_CO2[,"time"],"obsRes_CO2"] <- obsRes_CO2[,"Res_CO2"] 


