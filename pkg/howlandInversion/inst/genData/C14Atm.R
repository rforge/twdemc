# TODO: Add comment
# 
# Author: twutz
###############################################################################


delta14Catm <- read.csv(file.path("inst","genData","C14Atm.csv"))
str(delta14Catm)
names(delta14Catm) <- c("yr","delta14C")
plot(delta14C ~ yr, data=delta14Catm)

#plot(delta2iR14C(delta14C) ~ yr, data=delta14Catm) # atomic ratio as multiple of standard 1950  

save(delta14Catm,file=file.path("data","delta14Catm.RData"))
