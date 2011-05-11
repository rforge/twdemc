#install.packages("FME")
require(FME)
#?FME

windows(width=3.4,height=2.8, pointsize=10, record=TRUE)
par( las=1 )					#also y axis labels horizontal
par(mar=c(2.2,2.8,1,0)+0.3 )  #margins
par(tck=0.02 )				#axe-tick length inside plots             
par(mgp=c(1.4,0.2,0) )  #positioning of axis title, axis labels, axis


## =======================================================================
## ICBM
## =======================================================================
# assume we can mark litter inputs 
# and trace ratio of new Carbon in respiration and stocks

yini    <- c(
	Y = 1		##<< young carbon kg/ha
	,O=3		##<< old carbon kg/ha
	,YNew = 0
	,ONew = 0
)	 

pars    <- c(
	litter = 0.6	##<< kgC/ha/yr litter input
	,kY  = 1		##<< 1/yr decay rate young pool tvr of 1 year 
	,kO  = 1/40		##<< 1/yr decay rate of old pool
	,h = 0.15		##<< humification coefficient
)

yEq <- with(as.list(pars), litter/kY)
oEq <- with(as.list(pars), (h*yEq)/kO)

times   <- seq(0, 50, by = 1)	# years

# fluctuating inputs
inputs <- rnorm( length(times), pars["litter"], pars["litter"]*0.3 )
#fInput <- approxfun(times, inputs)
fInputConst <- function(times){ pars["litter"] }
plot( inputs ~ times )

MICBM2 <- function(Time, State, Pars, fInput) {
	with(as.list(c(State, Pars)), {
			decY  <- kY*Y
			decO  <- kO*O
			#dY    <- +fInput(Time) -decY
			dY    <- +0 -decY
			dO    <- +h*decY - decO
			
			decYNew <- kY*YNew
			decONew <- kO*ONew
			dYNew    <- +fInput(Time) -decYNew
			dONew    <- +h*decYNew - decONew
			
			resp  <- (1-h)*decY + decO
			respNew <- (1-h)*decYNew + decONew
			
			return(list(c(dY,dO,dYNew,dONew)
					,c(resp=resp,respNew=respNew)))
		})
}

# Solving the model
out     <- ode(func = MICBM2, y = yini, parms = pars, times = times, fInput=fInput)
tail(out)

# Plotting the results
matplot(out[,1], out[,2:5], type = "l", xlab = "time (yr)", ylab = "",	main = "ICBM", lwd = 2)
mtext("C (kg/ha/yr)",2,2.2, las=0)
legend("topright", c("Young","Old","Young New","Old New"), col = 1:4, lty = 1:4, inset=c(0.05,0.15))

# Plotting the results
matplot(out[,1], cbind( 
	out[,c("respNew")]/rowSums(out[,c("resp","respNew")])
	,rowSums(out[,c("YNew","ONew")])/rowSums(out[,c("Y","O","YNew","ONew")])
)
, type = "l", xlab = "time (yr)", ylab = "",	main = "ICBM", lwd = 2)

mtext("Proportion of new (labeled) C",2,2, las=0)
legend("bottomright", c("Respiration","Stocks"), col = 1:4, lty = 1:4, inset=c(0.05,0.15))
