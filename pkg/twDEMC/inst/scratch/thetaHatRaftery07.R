.tmp.estimateLMax.Raftery <- function(){
	# estimate lmax according to Raftery07: l_{max} = \bar{l} + Var{l}
	# does not work out, varL is much too high for strongly skewed distr.
	# need to thin because of autocorrlation before estimating the variance
	sscT <- squeeze( ssc, length.out=effectiveSize(mcl) )
	logDen <- stackChains(concatPops(sscT)$logDen)
	varL <- apply(logDen,2,var)
	meanL <- apply(logDen,2,mean)
	lMax1 <- meanL + varL
}
