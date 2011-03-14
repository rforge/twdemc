     \name{dsKinrespTut}
     \docType{data}
     \alias{dsKinrespTut}
     \title{Time series of soil respiration.}
     \description{}
     \usage{dsKinrespTut}
     \format{data.frame with columns \describe{
		\item{suite}{Grouping of experiments due to same experimental setup. }
		\item{experiment}{Identifier of experiment. }
		\item{replicate}{Replicate within experiment. }
		\item{resp}{Respiration in µg respired C/g Soil/hour}
		\item{time}{Time in hours}
		}
     }
     \source{subset from dataset of Sergey Blagodatsky}
     \references{
     	Wutzler et al.: Soil microbial biomass and its activity estimated by kinetic respiration analysis – Statistical guidelines., in prep.
     }
     \keyword{datasets}
     \examples{
     	data(dsKinrespTut)
     	plot( resp~time, col=replicate, data=subset(dsKinrespTut, experiment=="9"), xlab="Time (hours)", ylab="Respiration (µg Resp-C/g Soil/hour)" )
     }
     
