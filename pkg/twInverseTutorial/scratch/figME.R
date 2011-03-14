

windows(width=3.4,height=2.8, pointsize=10, record=TRUE)
par( las=1 )					#also y axis labels horizontal
par(mar=c(2.2,2.5,0,0)+0.3 )  #margins
par(tck=0.02 )				#axe-tick length inside plots             
par(mgp=c(1.4,0.2,0) )  #positioning of axis title, axis labels, axis

pHead<-pHead0<-0.7
n=30
y=0:n
#?pbinom
p<-dbinom(y,n,pHead)

plot(p~y,ylab="", xlab="Number of heads" )
mtext("Observation probability (pHead=0.7)",2,2,las=0)
library(ggplot2)
qplot(y,p,xlab="Number of heads",ylab="Observation probability (pHead=0.7)")

library(MASS)	#fitdistr
#?fitdistr
pHead <- seq(0.1,0.9,length.out=8*10+1)
y <- pHead0*n
p<-0.5-abs(pbinom(y,n,pHead)-0.5)
plot(p~pHead)
qplot(pHead,p,xlab="pHead",ylab="Parameter probability (nHead=21)",geom="line" )


