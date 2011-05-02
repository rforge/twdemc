## =======================================================================
## Example1: Predator-Prey Lotka-Volterra model
## =======================================================================

LVmod <- function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
			Ingestion    <- rIng  * Prey*Predator
			GrowthPrey   <- rGrow * Prey*(1-Prey/K)
			MortPredator <- rMort * Predator
			
			dPrey        <- GrowthPrey - Ingestion
			dPredator    <- Ingestion*assEff - MortPredator
			
			return(list(c(dPrey, dPredator)))
		})
}

pars    <- c(rIng   = 0.2,    # /day, rate of ingestion
	rGrow  = 1.0,    # /day, growth rate of prey
	rMort  = 0.2 ,   # /day, mortality rate of predator
	assEff = 0.5,    # -, assimilation efficiency
	K      = 10)     # mmol/m3, carrying capacity

yini    <- c(Prey = 1, Predator = 2)
times   <- seq(0, 200, by = 1)
out     <- ode(func = LVmod, y = yini, parms = pars, times = times)
summary(out)

matplot(out[,1], out[,2:3], type = "l", xlab = "time", ylab = "Conc",
	main = "Lotka-Volterra", lwd = 2)
legend("topright", c("prey", "predator"), col = 1:2, lty = 1:2)

## =======================================================================
## Example2: Substrate-Producer-Consumer Lotka-Volterra model
## =======================================================================

## Note:
## Function sigimp passed as an argument (input) to model
##   (see also lsoda and rk examples)

SPCmod <- function(t, x, parms, input)  {
	with(as.list(c(parms, x)), {
			import <- input(t)
			dS <- import - b*S*P + g*C    # substrate
			dP <- c*S*P  - d*C*P          # producer
			dC <- e*P*C  - f*C            # consumer
			res <- c(dS, dP, dC)
			list(res)
		})
}

## The parameters 
parms <- c(b = 0.001, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)

## vector of timesteps
times <- seq(0, 200, length = 101)

## external signal with rectangle impulse
signal <- as.data.frame(list(times = times,
		import = rep(0, length(times))))

signal$import[signal$times >= 10 & signal$times <= 11] <- 0.2

sigimp <- approxfun(signal$times, signal$import, rule = 2)


## Start values for steady state
xstart <- c(S = 1, P = 1, C = 1)

## Solve model
out <- ode(y = xstart,times = times,
	func = SPCmod, parms, input = sigimp)

## Default plot method
plot(out, type = "l")

## User specified plotting
mf <- par(mfrow = c(1, 2))
matplot(out[,1], out[,2:4], type = "l", xlab = "time", ylab = "state")
legend("topright", col = 1:3, lty = 1:3, legend = c("S", "P", "C"))
plot(out[,"P"], out[,"C"], type = "l", lwd = 2, xlab = "producer",
	ylab = "consumer")
par(mfrow = mf)
