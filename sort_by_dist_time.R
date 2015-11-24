#author: Lisa-Marie Harrison

#------------------------- sort by distance or time -------------------------------#
setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
glm.full <- read.csv(file = "fullCTD.csv", header= T) 
glm.full[which(glm.full$fluoro == -9), c(9:16)] <- NA
attach(glm.full)

#find absolute distance from station 2 for all stations
library(lmap)
stn.dist <- sqrt((lat - lat[1])^2 + (long - long[1])^2)
dist.from.stn2 <- rep(stn.dist, 1, each = 125)
glm.full <- cbind(glm.full, dist.from.stn2)

#find latitude and longitude distances from stn 2
lat.dist <- sqrt((lat - lat[1])^2)
lat.from.stn2 <- rep(lat.dist, 1, each = 125)
long.dist <- sqrt((long - long[1])^2)
long.from.stn2 <- rep(long.dist, 1, each = 125)
glm.full <- cbind(glm.full, lat.from.stn2, long.from.stn2)

#find time since start
ctd.info <- read.csv(file = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West raw data/Echoview/test/ctd_information.csv")
ctd.info <- ctd.info[-1, ] #rm station 1
time.from.start <- ctd.info$day - ctd.info$day[1]
for(i in 1:117){
  if (ctd.info$month[i] == 2) time.from.start[i] = 21 + ctd.info$day[i]
}
time.from.start <- rep(time.from.start, 1, each = 125)
glm.full <- cbind(glm.full, time.from.start)

#sort by distance from station or time since start of survey
sort.by.dist <- glm.full[order(glm.full$dist.from.stn2), ]
sort.by.time <- glm.full
sort.by.lat  <- glm.full[order(glm.full$lat.from.stn2, glm.full$long.from.stn2), ]
sort.by.long  <- glm.full[order(glm.full$long.from.stn2, glm.full$lat.from.stn2), ]


plot(sort.by.dist$l.fluoro[sort.by.dist$profile.depth == 40])
or <- 0
for (i in unique(glm.full$profile.depth)){
  a <- na.omit(sort.by.dist$l.fluoro[sort.by.dist$profile.depth == i])
  or[i] <- ar(a)$order
}
or <- na.omit(or)

ord <- 0
for (i in unique(glm.full$stn)){
  a <- na.omit(glm.full$l.fluoro[glm.full$stn == i])
  ord[i] <- ar(a)$order
}






