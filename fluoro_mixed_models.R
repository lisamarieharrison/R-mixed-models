setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.full <- read.table(file = "glm.full.txt", header= T)
stn.info <- read.table(file = "station_information.txt", skip = 1, fill = T)
attach(glm.full)
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

#------------------------ MIXED MODEL WITH ALL VARIABLES - NO CORRELATION TERMS --------------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + profile.depth + bot_nit + sal + temp + bot_oxy 
                    + bot_sil + ice.free.days + bot_phos, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(bot_nit) + spl(sal) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days) + spl(bot_phos)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)


#plot residuals for a certain station
s <- 40
plot(as.numeric(glm.full.sort$profile.depth[glm.full.sort$stn == s]), 
     fluoro.mm$residuals[glm.full.sort$stn == s], xlab = "profile depth (m)", ylab = "residuals")
title(bquote("Plot of residuals at station" ~ .(s)~""))


#plot the residuals at 40m as a bubble plot by latitude and longitude
d <- 40
l.fluoro.bubble <- fluoro.mm$residuals[glm.full.sort$profile.depth == d]
rad   <- sqrt(abs(l.fluoro.bubble)/pi) #set size = radius rather than area of circle
lat   <- glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d]
long  <- glm.full.sort$stn.long.full[glm.full.sort$profile.depth == d]
point_col <- rep(0, length(l.fluoro.bubble))
point_col[l.fluoro.bubble <= 0] <- "#0000FF"
point_col[l.fluoro.bubble > 0] <- "#FF0000"
symbols(x = long, y = lat, circles = rad, inches = 0.2, fg = point_col, lwd = 2.5)
title(bquote("Residuals by latitude/longitude at " ~ .(d)~"m depth"))
legend("topleft", c("negative", "positive"), col = c("blue", "red"), lwd = 2.5, bty = "n", cex = 0.75)


#correlogram to assess spatial auto-correlation
d <- 40
correlog <- spline.correlog(x = glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d], 
                            y = glm.full.sort$stn.long.full[glm.full.sort$profile.depth == d], 
                            z = fluoro.mm$residuals[glm.full.sort$profile.depth == d], na.rm = T)
plot.spline.correlog(correlog)
title(bquote("Correlogram at " ~ .(d)~"m depth"))


#variogram at depth 40
d <- 40
v1 <- variog(coords = cbind(glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d], 
                            glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d]), 
             data = fluoro.mm$residuals[glm.full.sort$profile.depth == d])
plot(v1, main = bquote("Semivariogram at " ~ .(d)~"m depth"), pch = 19)


#------------------- MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES - NO CORRELATION TERMS --------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#-------------- MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL DEPTH CORRELATION -------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days), rcov =~ stn:profile.depth
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)


#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL AR(1) DEPTH CORRELATION ------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days), rcov =~ stn:ar1(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL AR(2) DEPTH CORRELATION ------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days), rcov =~ stn:ar2(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL MA DEPTH CORRELATION ------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + bot_oxy + bot_sil + ice.free.days
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) 
                    , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:ma(profile.depth), 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)


#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL EXP DEPTH CORRELATION ------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + bot_oxy + bot_sil + ice.free.days
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) 
                    , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:exp(profile.depth), 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)


#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RANDOM MA DEPTH TERM ------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + bot_oxy + bot_sil + ice.free.days
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) + stn:ma(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#------------ MIXED MODEL WITH RANDOM AR(2) DEPTH TERM ------------#

#make a cut down data set with only every 4th station included 
glm.full.cut <- glm.full.sort[1, ]
for(i in 1:nrow(glm.full.sort)){
  if(any(seq(2, 118, , by = 4) == as.numeric(glm.full.sort$stn[i])))
    glm.full.cut <- rbind(glm.full.cut, glm.full.sort[i, ])
}
glm.full.cut <- glm.full.cut[-1, ]
glm.full.cut$profile.depth <- as.factor(glm.full.cut$profile.depth)


fluoro.mm <- asreml(fixed = l.fluoro ~ + par + bot_nit + temp + bot_oxy + bot_sil + ice.free.days
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(bot_oxy) + spl(bot_sil) + spl(par)
                    + spl(temp) + spl(ice.free.days) + stn:ar2(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.cut, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)



#plot residuals for a certain station
s <- 10
plot(as.numeric(levels(glm.full.cut$profile.depth[glm.full.cut$stn == s])), 
     fluoro.mm$residuals[glm.full.cut$stn == s], xlab = "profile depth (m)", ylab = "residuals")
title(bquote("Plot of residuals at station" ~ .(s)~""))


#plot the residuals at 40m as a bubble plot by latitude and longitude
d <- 38
l.fluoro.bubble <- fluoro.mm$residuals[glm.full.cut$profile.depth == d]
rad   <- sqrt(abs(l.fluoro.bubble)/pi) #set size = radius rather than area of circle
lat   <- glm.full.cut$stn.lat.full[glm.full.cut$profile.depth == d]
long  <- glm.full.cut$stn.long.full[glm.full.cut$profile.depth == d]
point_col <- rep(0, length(l.fluoro.bubble))
point_col[l.fluoro.bubble <= 0] <- "#0000FF"
point_col[l.fluoro.bubble > 0] <- "#FF0000"
symbols(x = long, y = lat, circles = rad, inches = 0.2, fg = point_col, lwd = 2.5)
title(bquote("Residuals by latitude/longitude at " ~ .(d)~"m depth"))
legend("topleft", c("negative", "positive"), col = c("blue", "red"), lwd = 2.5, bty = "n", cex = 0.75)


#correlogram to assess spatial auto-correlation
d <- 40
correlog <- spline.correlog(x = glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d], 
                            y = glm.full.sort$stn.long.full[glm.full.sort$profile.depth == d], 
                            z = fluoro.mm$residuals[glm.full.sort$profile.depth == d], na.rm = T)
plot.spline.correlog(correlog)
title(bquote("Correlogram at " ~ .(d)~"m depth"))


#variogram at depth 40
d <- 40
v1 <- variog(coords = cbind(glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d], 
                            glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d]), 
             data = fluoro.mm$residuals[glm.full.sort$profile.depth == d])
plot(v1, main = bquote("Semivariogram at " ~ .(d)~"m depth"), pch = 19)





#------------------ AUTOCORRELATION OF RESIDUALS FOR MODEL WITH NO AR TERM ---------------------#


fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + bot_oxy + bot_sil + ice.free.days
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) 
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)


#autocorrelation function for residuals at a station
s <- 100
a <- acf(fluoro.mm$residuals[glm.full.sort$stn == s], na.action = na.pass, lag.max = 125, main = "")
title(bquote("ACF function on residuals at station" ~ .(s)~""))



#------------------- MIXED MODEL RUN WITHOUT TRANSFORMATION TO LOG FLUORO ----------------------#


fluoro.mm <- asreml(fixed = fluoro ~ par + profile.depth + bot_nit + sal + temp + bot_oxy 
                    + bot_sil + ice.free.days + spl(bot_phos), random = ~ stn 
                    + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) +spl(sal) + spl(bot_phos)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)



#plot residuals for a certain station
s <- 40
plot(as.numeric(glm.full.sort$profile.depth[glm.full.sort$stn == s]), 
     fluoro.mm$residuals[glm.full.sort$stn == s], xlab = "profile depth (m)", ylab = "residuals")
title(bquote("Plot of residuals at station" ~ .(s)~""))


#plot the residuals at 40m as a bubble plot by latitude and longitude
d <- 40
l.fluoro.bubble <- fluoro.mm$residuals[glm.full.sort$profile.depth == d]
rad   <- sqrt(abs(l.fluoro.bubble)/pi) #set size = radius rather than area of circle
lat   <- glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d]
long  <- glm.full.sort$stn.long.full[glm.full.sort$profile.depth == d]
point_col <- rep(0, length(l.fluoro.bubble))
point_col[l.fluoro.bubble <= 0] <- "#0000FF"
point_col[l.fluoro.bubble > 0] <- "#FF0000"
symbols(x = long, y = lat, circles = rad, inches = 0.2, fg = point_col, lwd = 2.5)
title(bquote("Residuals by latitude/longitude at " ~ .(d)~"m depth"))
legend("topleft", c("negative", "positive"), col = c("blue", "red"), lwd = 2.5, bty = "n", cex = 0.75)


#correlogram to assess spatial auto-correlation
d <- 40
correlog <- spline.correlog(x = glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d], 
                            y = glm.full.sort$stn.long.full[glm.full.sort$profile.depth == d], 
                            z = fluoro.mm$residuals[glm.full.sort$profile.depth == d], na.rm = T)
plot.spline.correlog(correlog)
title(bquote("Correlogram at " ~ .(d)~"m depth"))


#variogram at depth 40
d <- 40
v1 <- variog(coords = cbind(glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d], 
                            glm.full.sort$stn.lat.full[glm.full.sort$profile.depth == d]), 
             data = fluoro.mm$residuals[glm.full.sort$profile.depth == d])
plot(v1, main = bquote("Semivariogram at " ~ .(d)~"m depth"), pch = 19)




#--------------------- MIXED MODEL WITH SQUARE ROOT TRANSFORMATION -------------------------#


fluoro.mm <- asreml(fixed = sqrt(fluoro) ~ par + profile.depth + bot_nit + sal + temp + bot_oxy 
                    + bot_sil + ice.free.days + spl(bot_phos), random = ~ stn 
                    + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) +spl(sal) + spl(bot_phos)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)


#--------------------- MIXED MODEL WITH GAUSSIAN DEPTH ERROR STRUCTURE -------------------------#


fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + sal + temp + bot_oxy + bot_sil + ice.free.days
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) +spl(sal), rcov =~ stn:gau(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)


#-------------------------------- MIXED MODEL WITH LATITUDE/LONGITUDE RANDOM FACTORS ------------------------------------#

fluoro.mm <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn + stn.lat.full + stn.long.full
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)




