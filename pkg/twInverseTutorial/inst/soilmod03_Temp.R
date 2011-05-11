## =======================================================================
## ICBM
## =======================================================================

tempFac <- function(temp, Q10=2, temp0=10){
	Q10^((temp-temp0)/10)	
}
times   <- seq(0, 200, by = 1)	# years
temp <- seq(10,length.out=length(times),by=5/100)
tempNoise <- temp + rnorm(length(temp), sd=2)
fTemp <- approxfun(times, tempNoise)

plot( tempFac(fTemp(times)) ~ times )

# incorporate temperature dependence in ICBM model


	 

pars    <- c(
	litter = 0.6	##<< kgC/ha/yr litter input
	,kY  = 1		##<< 1/yr decay rate young pool tvr of 1 year at temp0
	,kO  = 1/40		##<< 1/yr decay rate of old pool at temp0
	,h = 0.15		##<< humification coefficient
)

yEq <- with(as.list(pars), litter/kY)
oEq <- with(as.list(pars), (h*yEq)/kO)

yini    <- c(
	Y = yEq		##<< young carbon kg/ha
	,O=oEq		##<< old carbon kg/ha
)

# fluctuating inputs
inputs <- rnorm( length(times), pars["litter"], pars["litter"]*0.3 )
fInput <- approxfun(times, inputs)
#plot( inputs ~ times )

MICBM3 <- function(Time, State, Pars, fInput, fTemp) {
	with(as.list(c(State, Pars)), {
			tFac <- tempFac(fTemp(Time))
			decY  <- tFac*kY*Y
			decO  <- tFac*kO*O
			#dY    <- +fInput(Time) -decY
			dY    <- +litter -decY
			dO    <- +h*decY - decO
			resp  <- (1-h)*decY + decO
			return(list(c(dY,dO),c(decY=decY,decO=decO,resp=resp,tFac=tFac)))
		})
}

# Solving the model
out     <- ode(func = MICBM3, y = yini, parms = pars, times = times, fInput=fInput, fTemp=fTemp)
tail(out)

# Plotting the results
matplot(out[,1], out[,2:3], type = "l", xlab = "time (yr)", ylab = "",	main = "ICBM", lwd = 2)
mtext("C (kg/ha/yr)",2,2.2, las=0)
legend("topright", c("Young","Old"), col = 1:2, lty = 1:2, inset=c(0.05,0.05))


