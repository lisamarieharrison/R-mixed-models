#runs additive mixed models in R using mgcv'
setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
glm.full <- read.csv(file = "fullCTD.csv", header= T) 
ctd <- read.csv(file = "distToIce.csv", header = T)
glm.full[which(glm.full$fluoro == -9), c(9:16)] <- NA
attach(glm.full)


#plot of fluoro by station
par(mfrow = c(1, 2))
glm.full$stn <- as.numeric(as.factor(glm.full$stn))
plot(glm.full$fluoro, -glm.full$profile.depth, col = glm.full$stn, xlab = "fluoro", ylab = "depth")
title("Fluoro")
plot(glm.full$l.fluoro, -glm.full$profile.depth, col = glm.full$stn, xlab = "l.fluoro", ylab = "depth")
title("log(fluoro)")

#plot of station numbers by location
plot(long, lat, col = "white", xlab = "longitude",
     ylab = "latitude", main = "station number by latitude/longitude")
text(long, lat, labels = c(2:118))


#plot of fluoro by station with station numbers 
glm.full$stn <- as.numeric(as.factor(glm.full$stn))
plot(glm.full$fluoro, -glm.full$profile.depth, col = "white", xlab = "depth", ylab = "fluoro")
title("Fluoro")
text(glm.full$fluoro, -glm.full$profile.depth, labels = glm.full$stn)


#plot of station numbers by location with high fluoro locations in red
par(mfrow = c(1, 1))
glm.full$stn[which(glm.full$fluoro > 5)]
plot(long, lat, col = "white", xlab = "longitude",
     ylab = "latitude", main = "station number by latitude/longitude")
text(long, lat, labels = c(2:118))
rd.stn <- unique(glm.full$stn[which(glm.full$fluoro > 5)])
text(long[rd.stn], lat[rd.stn], labels = rd.stn, col = "red")


#plots of fluoro against each variable
par(mfrow = c(2, 2))
plot(glm.full$profile.depth, glm.full$fluoro, xlab = "profile.depth", ylab = "fluoro", main = "Depth")
plot(glm.full$temp, glm.full$fluoro, xlab = "temp", ylab = "fluoro", main = "Temperature")
plot(glm.full$par, glm.full$fluoro, xlab = "par", ylab = "fluoro", main = "PAR")
plot(glm.full$oxy, glm.full$fluoro, xlab = "oxy", ylab = "fluoro", main = "Oxygen")


#plot variables for each fluoro group
par(mfrow = c(2, 2))
group2 <- unique(glm.full$stn[which(glm.full$fluoro >= 6)])
plot(glm.full$profile.depth, glm.full$fluoro, xlab = "profile.depth", ylab = "fluoro", main = "Depth")
points(glm.full$profile.depth[stn %in% group2], glm.full$fluoro[stn %in% group2], col = "red")
plot(glm.full$temp, glm.full$fluoro, xlab = "temp", ylab = "fluoro", main = "Temperature")
points(glm.full$temp[stn %in% group2], glm.full$fluoro[stn %in% group2], col = "red")
plot(glm.full$par, glm.full$fluoro, xlab = "par", ylab = "fluoro", main = "PAR")
points(glm.full$par[stn %in% group2], glm.full$fluoro[stn %in% group2], col = "red")
plot(glm.full$oxy, glm.full$fluoro, xlab = "oxy", ylab = "fluoro", main = "Oxygen")
points(glm.full$oxy[stn %in% group2], glm.full$fluoro[stn %in% group2], col = "red")

#plot profile depth against par with fluoro as colour
z <- na.omit(glm.full$fluoro)
zlen <- length(z)
levels <- seq(min(z),max(z),length.out = zlen)
col <- colorRampPalette(c("black","yellow"))(zlen)[rank(z)]
plot(glm.full$par, -glm.full$profile.depth, col = col, pch = 19)


#plot distance from ice edge by latitude/longitude
plot(ctd$longitude, ctd$latitude, col = "white", xlab = "longitude",
     ylab = "latitude", main = "distance to ice by latitude/longitude")
text(ctd$longitude, ctd$latitude, labels = round(ctd$distToIce), 0)

#plot fluoro against distance to ice edge 
plot(rep(ctd$distToIce, 1, each = 125), glm.full$fluoro[glm.full$stn != 1], xlab = "Distance to Ice", ylab = "fluoro")
title("fluoro against distance to ice edge")

#plot fluoro against depth with distance from ice edge <100km in red
plot(glm.full$profile.depth, glm.full$fluoro, xlab = "profile depth", ylab = "fluoro")
rd.stn <- unique(ctd$station[which(ctd$distToIce <100)])
points(glm.full$profile.depth[glm.full$stn %in% rd.stn], glm.full$fluoro[glm.full$stn %in% rd.stn], col = "red")

#plot distance from ice edge against par with fluoro > 6 in red
plot(rep(ctd$distToIce, 1, each = 125)[glm.full$profile.depth < 50], glm.full$par[glm.full$profile.depth < 50])
rd.stn <- unique(glm.full$stn[which(glm.full$fluoro > 6)])
points(rep(ctd$distToIce, 1, each = 125)[glm.full$profile.depth < 50 & glm.full$stn %in% rd.stn], glm.full$par[glm.full$profile.depth < 50 & glm.full$stn %in% rd.stn], col = "red")













