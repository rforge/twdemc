levels( ceCoef$compart )
iComp="stem"
res <- t(do.call(cbind, lapply( c("stem","branches","leaves","root"), function(iComp){
		print(iComp)
		iSpec="pine"
		iSpecs=
		tmp <- rbind( iComp, do.call( cbind, lapply( levels(ceCoef$species), function(iSpec){
						print(iSpec)
						sis = unique(ceCoef$si[ ceCoef$species==iSpec & ceCoef$compartment==iComp ])
						iSi = sis[1]
						rbind(iSpec, sapply( sis, function(iSi){
									print(iSi)
									c(levels(ceCoef$si)[iSi], format( calcCE(iSpec, age, iComp, iSi ), digits=2 ))
								}))
					})))
		tmp
	}) ))
colnames(res) <- c("Compartment","Species","Site index",paste(age,"yr",sep=""))
copy2clip(res, row.names=FALSE, quote=FALSE)

