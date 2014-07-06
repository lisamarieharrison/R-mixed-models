#runs mixed model on thinned data set using CTD data only

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
#usually use thinCTD.csv, fullCTD.csv is the full data set using CTD data only and padding missing depths with NA
glm.full <- read.csv(file = "rstnCTD.csv", header= T) 
library(asreml)
library(rgl)
library(akima)
library(gstat)
attach(glm.full)


#basic mixed model
fluoro.mm1 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + oxy, random = ~ stn  
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(oxy)  
                     , na.method.X = " include ", na.method.Y = "include", 
                     data = glm.full, maxiter = 50, workspace = 400000000)
summary(fluoro.mm1)
wald(fluoro.mm1)
plot(fluoro.mm1)

#3D plot of residuals for a single depth to assess latitude, longitude spatial autocorrelation
zgamma <- resid(fluoro.mm1)
mat <- cbind(abs(glm.full$lat[glm.full$profile.depth == 12]), glm.full$long[glm.full$profile.depth == 12], zgamma[glm.full$profile.depth == 12])
colnames(mat) <- c("lat", "long", "zgamma")
mat <- na.omit(mat)
s <- interp(mat[, 1], mat[, 2], mat[, 3])
s$z[s$z == Inf] <- NA

zlim <- range(s$z, na.rm = TRUE)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- terrain.colors(zlen,alpha=0) # height color lookup table
col <- colorlut[ s$z-zlim[1]+1 ] # assign colors to heights for each point

surface3d(s$x, s$y, s$z, color = col, back = "lines")
axes3d('x--', labels = round(s$x, 2), color = "black")
axes3d('y--', labels = round(s$y, 2), color = "black")
axes3d('z--', labels = round(range(s$z, na.rm = TRUE), 2), color = "black")

#semivariogram of residuals by depth for a station
s <- 2 #choose station number
dat <- na.omit(cbind(resid(fluoro.mm1)[glm.full$stn == s], glm.full$profile.depth[glm.full$stn == s]))
v1 <- variogram(list(dat[, 1]), locations = list(dat[, 2]))
plot(v1$dist, v1$gamma, pch = 19, xlab = "distance (depths)", type = "l", ylab = "semivariance")
title("Semivariogram for each station")

#compute a semivariogram by depth for each station and add to the plot
for(i in unique(stn)){
  dat <- na.omit(cbind(resid(fluoro.mm1)[glm.full$stn == i], glm.full$profile.depth[glm.full$stn == i]))
  assign(paste("v", i, sep = ""), variogram(list(dat[, 1]), locations = list(dat[, 2])))
  v1 <- variogram(list(dat[, 1]), locations = list(dat[, 2]))
  points(v1$dist, v1$gamma, pch = 19, type = "l")
}


#mixed model with simple correlation structure
fluoro.mm1 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + oxy, random = ~ stn  
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(oxy), rcov=~ stn:profile.depth  
                     , na.method.X = " include ", na.method.Y = "include", 
                     data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm1)
wald(fluoro.mm1)
plot(fluoro.mm1)


#mixed model with moving average correlation structure
fluoro.mm1 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + oxy, random = ~ stn  
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(oxy), rcov=~ stn:ma(profile.depth)
                     , na.method.X = " include ", na.method.Y = "include", 
                     data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm1)
wald(fluoro.mm1)
plot(fluoro.mm1)



#mixed model with correlation as random effect
depth.ar <- corGaus(form =~ stn|profile.depth)
fluoro.mm1 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + oxy, random = ~ stn  
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(oxy)
                     , na.method.X = " include ", na.method.Y = "include", rov =~ corSpher(form =~ profile.depth|stn),
                     data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm1)
wald(fluoro.mm1)
plot(fluoro.mm1)

#create depth factor for use in correlation matrices
profile.depth.f <- as.factor(profile.depth)
glm.full <- cbind(glm.full, profile.depth.f)

#basic mixed model with depth as a random effect
fluoro.mm1 <- asreml(fixed = l.fluoro ~ par + temp + oxy + profile.depth, random = ~ stn
                     + spl(par) + spl(temp) + spl(oxy) + spl(profile.depth)
                     , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:fa(profile.depth.f),
                     data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm1)
wald(fluoro.mm1)
plot(fluoro.mm1)                       


#asreml model with blocked splines
fluoro.mm1 <- asreml(fixed = l.fluoro ~ par + temp + profile.depth + oxy,
                    random =~ stn:spl(par) +  spl(par) + spl(profile.depth) + spl(temp) + spl(oxy),
                    rcov =~ id(units), data = glm.full, na.method.Y = "include", na.method.X = "include", maxiter = 30, workspace = 400000000)










