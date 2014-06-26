setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.full <- read.table(file = "glm.full.txt", header= T)
stn.info <- read.table(file = "station_information.txt", skip = 1, fill = T)
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

glm.full          <- cbind(glm.full, stn.lat.full, stn.long.full)
glm.full.sort     <- glm.full[order(stn, profile.depth), ]
glm.full.sort$stn <- as.factor(glm.full.sort$stn)
glm.full.cut$stn  <- as.factor(glm.full.cut$stn)

#remove observations from first test station 
glm.full.sort <- glm.full.sort[which(glm.full.sort$stn.long.full < 85), ]

glm.full.cut <- glm.full.sort[1, ]
for(i in 1:nrow(glm.full.sort)){
  if(any(seq(2, 118, , by = 4) == as.numeric(glm.full.sort$stn[i])))
  glm.full.cut <- rbind(glm.full.cut, glm.full.sort[i, ])
}
glm.full.cut <- glm.full.cut[-1, ]
glm.full.cut$profile.depth <- as.factor(glm.full.cut$profile.depth)

#----------------------------- STANDARD MIXED MODEL -------------------------------#

#all variables
fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(profile.depth) + spl(bot_nit) + spl(par) 
                    + spl(temp) + stn:ar1(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.cut, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#Drop the least significant variable one at a time, with the best model being:
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + sal, random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(profile.depth) + spl(par) 
                    + spl(sal), na.method.X = " include ", na.method.Y = "include", rcov =~,
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

plot(glm.full.sort$profile.depth[glm.full.sort$stn == 10], fluoro.mm$residuals[glm.full.sort$stn == 10])

#plot the residuals at 40m as a bubble plot by latitude and longitude
d <- 40
l.fluoro.bubble <- fluoro.mm$residuals[stn.long.full < 85 & profile.depth == d]
rad   <- sqrt(abs(resid)/pi) #set size = radius rather than area of circle
lat   <- stn.lat.full[profile.depth == d & stn.long.full < 85]
long  <- stn.long.full[profile.depth == d & stn.long.full < 85]
symbols(x = long, y = lat, circles = rad, inches = 0.3)
title("Residuals for depth = 40m")


#ar process by depth

ar.stn2 <- ar(na.omit(glm.full.sort$l.fluoro[glm.full.sort$stn == 2]))

a <- 0.1105
b <- 0.6277158
d <- 0.1979

ar.i <- rep(0, 125)
ar.i[2] <- l.fluoro[stn == 2 & profile.depth == 4]
for(i in 4:125){
  ar.i[i] <- a*ar.i[i - 1] + b*ar.i[i - 2] + d*ar.i[i - 3]
}

plot(l.fluoro[stn == 100])
points(resid[stn == 100], col = "green")

a <- matrix(0, ncol = 118, nrow = 125)
for(i in 2:118){
  a[, i] <- runmed(l.fluoro[stn == i], k = 5, endrule = "keep")
}
b <- c(as.vector(a[, 2:118]))


#Drop the least significant variable one at a time, with the best model being:
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + temp, random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(profile.depth) + spl(par) + spl(temp)
                    , na.method.X = " include ", na.method.Y = "include"),
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

plot(glm.full.sort$profile.depth[glm.full.sort$stn == 100], fluoro.mm$residuals[glm.full.sort$stn == 100])



#variogram for epths
library(gstat)
v1 <- variogram(profile.depth, 
             data = na.omit(l.fluoro))
plot(v1, main = "Variogram of all stations at 40m depth")

correlog(profile.depth, l.fluoro)
