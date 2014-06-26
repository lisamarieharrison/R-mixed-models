setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.full <- read.table(file = "glm.full.txt", header= T)
stn.info <- read.table(file = "station_information.txt", skip = 1, fill = T)
glm.full.cut <- read.table(file = "glm.full.cut.txt", header = T)
attach(glm.full)
library(asreml)
library(calibrate)
library(ncf)
library(nlme)

# --------------- ADDING A SINGLE LATITUDE/LONGITUDE FOR EACH STATION ---------------#

#single lat/long for each station from station information pdf
stn.lat  <- -(stn.info[, 15] + stn.info[, 16]/60) #convert from dms to decimal degrees
stn.long <- stn.info[, 18] + stn.info[, 19]/60  #convert from dms to decimal degrees
plot(stn.long[stn.long < 85], stn.lat[stn.long < 85], col = "white", xlab = "longitude",
     ylab = "latitude", main = "station number by latitude/longitude")
text(stn.long[stn.long < 85], stn.lat[stn.long < 85], labels = c(1:118))

#add noise around duplicated stations to create 118 unique lat/long coordinates
a <- which(duplicated(stn.long) == T) 
stn.long[a] <- stn.long[a] + 0.0000000001
b <- which(duplicated(stn.lat) == T) 
stn.lat[b] <- stn.lat[b] + 0.0000000001

stn.long.full <- 0
stn.lat.full <- 0
for(i in 1:14750){
  stn.long.full[i] <- stn.long[stn[i]]
  stn.lat.full[i] <- stn.lat[stn[i]]
}

glm.full <- cbind(glm.full, stn.lat.full, stn.long.full)
glm.full.sort <- glm.full[order(stn, profile.depth), ]
glm.full.sort$stn <- as.factor(glm.full.sort$stn)
glm.full.sort$profile.depth <- as.factor(glm.full.sort$profile.depth)
glm.full.cut$profile.depth  <- as.factor(glm.full.cut$profile.depth)

#remove observations from first test station 
glm.full.sort <- glm.full.sort[which(glm.full.sort$stn.long.full < 85), ]


#----------------------------- STANDARD MIXED MODEL -------------------------------#

##DO NOT CHANGE!!!!

#mixed model
fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + sal
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(sal)
                    + spl(temp) + stn:ar2(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#plot residuals by depth
plot(rep(seq(2, 250, by = 2), 117), fluoro.mm$residuals)

#plot residuals for a certain depth
d <- 100
plot(as.numeric(levels(glm.full.sort$profile.depth[glm.full.sort$stn == d])), 
     fluoro.mm$residuals[glm.full.sort$stn == d])

#plot the residuals at 40m as a bubble plot by latitude and longitude
d <- 100
l.fluoro.bubble <- fluoro.mm$residuals[glm.full.sort$profile.depth == d]
rad   <- sqrt(abs(l.fluoro.bubble)/pi) #set size = radius rather than area of circle
lat   <- glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d]
long  <- glm.full.sort$stn.long.full[glm.full.sort$profile.depth == d]
symbols(x = long, y = lat, circles = rad, inches = 0.2)
title("Residuals for depth = 40m")

#------------------------- MODEL SELECTION ON CUT DOWN MATRIX -----------------------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + bot_oxy
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(bot_oxy)
                    + spl(temp) + stn:ar2(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
fluoro.mm <- update(fluoro.mm)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)



#plot residuals for a certain depth
d <- 38
plot(as.numeric(levels(glm.full.cut$profile.depth[glm.full.cut$stn == d])), 
     fluoro.mm$residuals[glm.full.cut$stn == d])

#plot the residuals at 40m as a bubble plot by latitude and longitude
d <- 100
l.fluoro.bubble <- fluoro.mm$residuals[glm.full.cut$profile.depth == d]
rad   <- sqrt(abs(l.fluoro.bubble)/pi) #set size = radius rather than area of circle
lat   <- glm.full.cut$stn.lat.full[glm.full.cut$profile.depth == d]
long  <- glm.full.cut$stn.long.full[glm.full.cut$profile.depth == d]
symbols(x = long, y = lat, circles = rad, inches = 0.2)
title("Residuals for depth = 40m")






