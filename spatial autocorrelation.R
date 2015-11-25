#author: Lisa-Marie Harrison

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.full <- read.table(file = "glm.full.txt", header= T)
stn.info <- read.table(file = "station_information.txt", skip = 1, fill = T)
library(asreml)
library(calibrate)
library(ncf)
library(nlme)
library(geoR)

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

stn.lat.full  <- sort(rep(stn.lat, 125))
stn.long.full <- sort(rep(stn.long, 125))

glm.full <- cbind(glm.full, stn.lat.full, stn.long.full)
glm.full.sort <- glm.full[order(glm.full$profile.depth), ]


#----------------------------- STANDARD MIXED MODEL -------------------------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + temp + par , random = ~ 
                      + spl(profile.depth) + stn + spl(temp) + spl(par)  
                    + stn.lat.full + stn.long.full, data = glm.full.sort,
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#remove latitude and longitude
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + temp + par , random = ~ 
                      + spl(profile.depth) + stn + spl(temp) + spl(par), data = glm.full.sort,
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#add all variables
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + sal + temp + bot_oxy + 
                      bot_sil + bot_phos + bot_nit, random = ~ stn + spl(profile.depth)
                    + spl(bot_nit) + spl(par) + spl(sal) + spl(temp) + spl(bot_oxy)
                    + spl(bot_sil) + spl(bot_phos), na.method.X = "include", 
                    na.method.Y = "include", data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#drop bottled oxygen
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + sal + temp + 
                      bot_sil + bot_phos + bot_nit, random = ~ stn + spl(profile.depth)
                    + spl(bot_nit) + spl(par) + spl(sal) + spl(temp) 
                    + spl(bot_sil) + spl(bot_phos), na.method.X = "include", 
                    na.method.Y = "include", data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#drop bottled silicon
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + sal + temp 
                    + bot_phos + bot_nit, random = ~ stn + spl(profile.depth)
                    + spl(bot_nit) + spl(par) + spl(sal) + spl(temp) 
                    + spl(bot_phos), na.method.X = "include", 
                    na.method.Y = "include", data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#drop salinity
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + temp 
                    + bot_phos + bot_nit, random = ~ stn + spl(profile.depth)
                    + spl(bot_nit) + spl(par) + spl(temp) 
                    + spl(bot_phos), na.method.X = "include", 
                    na.method.Y = "include", data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#drop temperature
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par 
                    + bot_phos + bot_nit, random = ~ stn + spl(profile.depth)
                    + spl(bot_nit) + spl(par)
                    + spl(bot_phos), na.method.X = "include", 
                    na.method.Y = "include", data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#drop bottled nitrogen
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par 
                    + bot_phos, random = ~ stn + spl(profile.depth) + spl(par)
                    + spl(bot_phos), na.method.X = "include", 
                    na.method.Y = "include", data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

# --------------------------- ADD AUTOCORRELATION BY DEPTH ---------------------------#

#moving average autocorrelation term for depth
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par 
                    + bot_phos, random = ~ stn + spl(profile.depth) + spl(par)
                    + spl(bot_phos), na.method.X = "include", na.method.Y = "include", 
                    data = glm.full, rcov =~ stn:ma(profile.depth), maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#find mean l.fluoro at each depth to use as autocorrelation term
depth.ar <- 0
for(i in 1:length(unique(profile.depth))){
  depth.ar[i] <- mean(na.omit(l.fluoro[profile.depth == profile.depth[i]]))
}

plot(profile.depth, l.fluoro, pch = 19)
points(seq(1, 250, 2), depth.ar, col = "green", pch = 19)
title("Log fluoro by depth with average l.fluoro at each depth in green")

depth.ar.full <- rep(depth.ar, 118)

glm.full <- cbind(glm.full, depth.ar.full)
glm.full.sort <- glm.full[order(stn, depth.ar.full), ]

#using mean l.fluoro by depth as autocorrelation term
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par 
                    + bot_phos, random = ~ stn + spl(profile.depth) + spl(par)
                    + spl(bot_phos), na.method.X = "include", na.method.Y = "include", 
                    data = glm.full.sort, rcov =~ stn:depth.ar.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#using depth as a correlation term
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par 
                    + bot_phos, random = ~ stn + spl(profile.depth) + spl(par)
                    + spl(bot_phos), na.method.X = "include", na.method.Y = "include", 
                    data = glm.full, rcov =~ stn:profile.depth, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#ordering by depth rather than station, but not including a correlation term
glm.full.sort <- glm.full[order(profile.depth, stn.lat.full), ]
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par 
                    + bot_phos, random = ~ stn + spl(profile.depth) + spl(par)
                    + spl(bot_phos), na.method.X = "include", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#---------------------------- ADD SPATIAL AUTOCORRELATION ----------------------------#

#autocorrelation by longitude
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + temp + par , random = ~ 
                      + spl(profile.depth) + stn + spl(temp) + spl(par)  
                    + stn.lat.full + stn.long.full, data = glm.full.sort, rcov =~ profile.depth:ma(stn.long.full),
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#check for autocorrelation of residuals by depth
plot(profile.depth[stn == 100], fluoro.mm$residuals[stn == 100], 
     xlab = "profile.depth", ylab = "residuals", main = "Residuals for station 100")

#check for autocorrelation of residuals by distance
plot(stn.lat.full[profile.depth == 40], fluoro.mm$residuals[profile.depth == 40], 
     xlab = "latitude", ylab = "residuals", main = "Residuals for depth == 40")
plot(stn.long.full[profile.depth == 40], fluoro.mm$residuals[profile.depth == 40], 
     xlab = "longitude", ylab = "residuals", main = "Residuals for depth == 40")

d <- 40
correlog <- spline.correlog(x = stn.lat.full[profile.depth == d], 
                            y = stn.long.full[profile.depth == d], 
                            z = fluoro.mm$residuals[profile.depth == d], na.rm = T)
plot.spline.correlog(correlog)
main_title <- bquote("Correlogram at " ~ .(d)~"m depth")
title(main_title)


#autocorrelation by latitude
glm.full.sort <- glm.full[order(profile.depth, stn.lat.full), ]

fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + temp + par , random = ~ 
                      + spl(profile.depth) + stn + spl(temp) + spl(par)  
                    + stn.lat.full + stn.long.full, data = glm.full.sort, rcov =~ profile.depth:ma(stn.lat.full),
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#check for autocorrelation of residuals by depth
plot(profile.depth[stn == 100], fluoro.mm$residuals[stn == 100], 
     xlab = "profile.depth", ylab = "residuals", main = "Residuals for station 100")

#check for autocorrelation of residuals by distance
plot(stn.lat.full[profile.depth == 40], fluoro.mm$residuals[profile.depth == 40], 
     xlab = "latitude", ylab = "residuals", main = "Residuals for depth == 40")
plot(stn.long.full[profile.depth == 40], fluoro.mm$residuals[profile.depth == 40], 
     xlab = "longitude", ylab = "residuals", main = "Residuals for depth == 40")

d <- 40
correlog <- spline.correlog(x = stn.lat.full[profile.depth == d], 
                            y = stn.long.full[profile.depth == d], 
                            z = fluoro.mm$residuals[profile.depth == d], na.rm = T)
plot.spline.correlog(correlog)
main_title <- bquote("Correlogram at " ~ .(d)~"m depth")
title(main_title)



#autocorrelation by latitude and longitude
glm.full$profile.depth <- as.factor(profile.depth)
glm.full$stn <- as.factor(stn)
glm.full.sort <- glm.full[order(profile.depth, stn.lat.full, stn.long.full), ]
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + temp + par , random = ~ 
                    + spl(profile.depth) + stn + spl(temp) + spl(par)  
                    + stn.lat.full + stn.long.full, data = glm.full, 
                    rcov =~ stn:ma(profile.depth),
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)


d <- 40
v1 <- variog(coords = cbind(stn.lat.full[profile.depth == d], stn.long.full[profile.depth == d]), 
             data = na.omit(l.fluoro[profile.depth == d]))
plot(v1, main = "Variogram of all stations at 40m depth")

variog.data <- cbind(stn.lat.full, stn.long.full, l.fluoro)
variog.data <- na.omit(variog.data)
v1 <- variog(coords = variog.data[, 1:2], 
             data = variog.data[, 3])
plot(v1)
