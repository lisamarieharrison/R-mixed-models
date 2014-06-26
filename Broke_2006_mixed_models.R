setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.full <- read.csv(file = "glm_full.csv", header= T)
glm.full.cut <- read.csv(file = "glm_full_cut.csv", header= T)
attach(glm.full)
library(asreml)
library(calibrate)
library(ncf)
library(nlme)
library(geoR)


#-------------------------------- EXPLORATORY PLOTS -------------------------------------------#

#plot of station number by latitude/longitude 
plot(unique(stn.long.full), unique(stn.lat.full), col = "white", xlab = "longitude",
     ylab = "latitude", main = "station number by latitude/longitude")
text(unique(stn.long.full), unique(stn.lat.full), labels = c(1:118))
range(profile.depth[stn == 1])

#plot of time since start of survey at each station
plot(stn.long[2:118], stn.lat[2:118], col = "white", xlab = "longitude",
     ylab = "latitude", main = "time since start of survey at each CTD station")
text(stn.long[2:118], stn.lat[2:118], labels = time[profile.depth == 2][2:118])


#plot l.fluoro against time for 40m depth
plot(time[profile.depth == 40], l.fluoro[profile.depth == 40], xlab = "days since start of survey",
     ylab = "l.fluoro", main = "log fluoro at 40m depth over the survey")


#plot l.fluoro against time for all depths
rbPal <- colorRampPalette(c('red','blue'))
dat_cols <- rbPal(125)[as.numeric(cut(profile.depth,breaks = 125))]
plot(c(time, c(51:60)), c(l.fluoro, rep(NA, 10)), col = dat_cols, ylab = "l.fluoro",
     xlab = "days since start of survey", main = "log fluro for entire survey", pch = 16)
legend("right", legend = seq(0, max(profile.depth), 20), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.6)


#plot l.fluoro against latitude/longitude for 40m depth - bubble plot
l.fluoro.bubble <- l.fluoro[stn.long.full < 85 & profile.depth == 40]
rad   <- sqrt(abs(resid)/pi) #set size = radius rather than area of circle
lat   <- stn.lat.full[profile.depth == 40 & stn.long.full < 85]
long  <- stn.long.full[profile.depth == 40 & stn.long.full < 85]
symbols(x = long, y = lat, circles = rad, inches = 0.3, xlab = "longitude", ylab = "latitude")
title("Log fluoro for depth = 40m")



#semivariograms to assess spatial autocorrelation
par(mfrow = c(2, 2))
par(oma = c(4, 4, 4, 0))
par(mar = c(2, 2, 1, 1))

d <- 110
v1 <- variog(coords = cbind(na.omit(glm.full)$stn.lat.full[na.omit(glm.full)$profile.depth == d], 
                            na.omit(glm.full)$stn.lat.full[na.omit(glm.full)$profile.depth == d]), 
             data = na.omit(glm.full)$l.fluoro[na.omit(glm.full)$profile.depth == d])
plot(v1, main = bquote("" ~ .(d)~"m depth"), pch = 19)

title("Semivariograms of log fluoro at various depths", outer = T)
mtext('distance', side = 1, outer = TRUE, line = 2)
mtext('semivariance', side = 2, outer = TRUE, line = 2)



#---------------------- PLOTS OF L.FLUORO AGAINST EXPLANATORY VARIABLES --------------------------#

#graphical parameters
par(xpd=T, mar=c(4,4,4,6))
rbPal <- colorRampPalette(c('red','blue'))
dat_cols <- rbPal(125)[as.numeric(cut(profile.depth,breaks = 125))]

#plot of bottled nitrogen against log fluoro
plot(bot_nit, l.fluoro, xlab = "bottled nitrogen", ylab = "l.fluoro", col = dat_cols, pch = 19)
title("log fluoro against bottled nitrogen")
legend(37, 4, legend = seq(0, max(profile.depth), 20), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.7, title = "depth (m)")


#plot of bottled silicon against log fluoro
plot(bot_sil, l.fluoro, xlab = "bottled silicon", ylab = "l.fluoro", col = dat_cols, pch = 19)
title("log fluoro against bottled silicon")
legend(103, 4, legend = seq(0, max(profile.depth), 20), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.7, title = "depth (m)")


#plot of bottled oxygen against log fluoro
plot(bot_oxy, l.fluoro, xlab = "bottled oxygen", ylab = "l.fluoro", col = dat_cols, pch = 19)
title("log fluoro against bottled oxygen")
legend(445, 4, legend = seq(0, max(profile.depth), 20), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.7, title = "depth (m)")

#plot of bottled phosphorus against log fluoro
plot(bot_phos, l.fluoro, xlab = "bottled phosphorus", ylab = "l.fluoro", col = dat_cols, pch = 19)
title("log fluoro against bottled phosphorus")
legend(2.5, 4, legend = seq(0, max(profile.depth), 20), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.7, title = "depth (m)")


#plot of temperature against log fluoro
plot(temp, l.fluoro, xlab = "temperature", ylab = "l.fluoro", col = dat_cols, pch = 19)
title("log fluoro against temperature")
legend(2, 4, legend = seq(0, max(profile.depth), 20), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.7, title = "depth (m)")


#plot of salinity against log fluoro
plot(sal, l.fluoro, xlab = "salinity", ylab = "l.fluoro", col = dat_cols, pch = 19)
title("log fluoro against salinity")
legend(35.1, 4, legend = seq(0, max(profile.depth), 20), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.7, title = "depth (m)")

#plot of ice free days against log fluoro
plot(ice.free.days[profile.depth == 40], l.fluoro[profile.depth == 40], xlab = "ice free days", ylab = "l.fluoro", pch = 19)
title("log fluoro against ice free days")

#plot of par against log fluoro
plot(par, l.fluoro, xlab = "par", ylab = "l.fluoro", col = dat_cols, pch = 19)
title("log fluoro against par")
legend(2, 4, legend = seq(0, max(profile.depth), 20), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.7, title = "depth (m)")

#plot of par against profile depth
par(xpd=F, mar=c(4,4,4,4))
plot(profile.depth, par, main = "profile depth against par")



#--------------------------------- 1. BASIC MIXED MODEL-----------------------------------#


fluoro.mm1 <- asreml(fixed = fluoro ~ par + profile.depth + temp + bot_oxy + bot_phos 
                     + bot_sil + ice.free.days, random = ~ stn 
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days) + spl(bot_phos)
                     , na.method.X = " include ", na.method.Y = "include", 
                     data = glm.full, maxiter = 30, workspace = 400000000)



#----------------------------- 2. SQUARE ROOT TRANSFORMATION --------------------------------#


fluoro.mm2 <- asreml(fixed = sqrt(fluoro) ~ par + profile.depth + bot_nit + sal + temp + bot_oxy 
                     + bot_sil + ice.free.days + bot_phos, random = ~ stn 
                     + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                     + spl(temp) + spl(ice.free.days) +spl(sal) + spl(bot_phos)
                     , na.method.X = " include ", na.method.Y = "include", 
                     data = glm.full, maxiter = 30, workspace = 400000000)


#----------------------------------- 3. LOG TRANSFORMATION ------------------------------------#

fluoro.mm3 <- asreml(fixed = l.fluoro ~ par + profile.depth + bot_nit + sal + temp + bot_oxy 
                    + bot_sil + ice.free.days + bot_phos, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(bot_nit) + spl(sal) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days) + spl(bot_phos)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full, maxiter = 30, workspace = 400000000)


#----------------------------------- 4. MODEL SELECTION ---------------------------------------#

fluoro.mm4 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full, maxiter = 30, workspace = 400000000)


#---------------====---- 5. LATITUDE/LONGITUDE RANDOM FACTORS ----------------------------------#

fluoro.mm5 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                     + bot_sil + ice.free.days, random = ~ stn + stn.lat.full + stn.long.full
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days)
                     , na.method.X = " include ", na.method.Y = "include", 
                     data = glm.full, maxiter = 30, workspace = 400000000)

#-------------------------- AUTOCORRELATION OF RESIDUALS --------------------------------------#
s <- 100
a <- acf(fluoro.mm$residuals[glm.full$stn == s], na.action = na.pass, lag.max = 125, main = "")
title(bquote("ACF function on residuals at station" ~ .(s)~""))



#------------------------- 6. RESIDUAL DEPTH ERROR STRUCTURE -----------------------------------#

fluoro.mm6 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days), rcov =~ stn:profile.depth
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full, maxiter = 30, workspace = 400000000)


#------------------------------ 7. RESIDUAL AR(1) DEPTH CORRELATION ----------------------------#

glm.full$profile.depth <- as.factor(glm.full$profile.depth)
fluoro.mm7 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                     + bot_sil + ice.free.days, random = ~ stn 
                     + spl(par) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days) 
                     , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:ar1(profile.depth),
                     data = glm.full, maxiter = 30, workspace = 400000000)



#------------------------------ 8. RESIDUAL AR(2) DEPTH CORRELATION ----------------------------#

glm.full$profile.depth <- as.factor(glm.full$profile.depth)
fluoro.mm8 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                     + bot_sil + ice.free.days, random = ~ stn 
                     + spl(par) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days) 
                     , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:ar2(profile.depth),
                     data = glm.full, maxiter = 30, workspace = 400000000)

#------------------------------ 9. RESIDUAL MA DEPTH CORRELATION --------------------------------#

glm.full$profile.depth <- as.numeric(glm.full$profile.depth)
fluoro.mm9 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                     + bot_sil + ice.free.days, random = ~ stn 
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days)
                     , na.method.X = " include ", na.method.Y = "include", rov =~ stn:ma(profile.depth), 
                     data = glm.full, maxiter = 30, workspace = 400000000)




#------------------------------ 10. RESIDUAL EXP DEPTH CORRELATION -----------------------------#

fluoro.mm10 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                     + bot_sil + ice.free.days, random = ~ stn 
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days)
                     , na.method.X = " include ", na.method.Y = "include", rcov =~ stn: exp(profile.depth), 
                     data = glm.full, maxiter = 30, workspace = 400000000)


#-------------------------------------- 11. RANDOM MA DEPTH TERM -------------------------------#

fluoro.mm11 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                     + bot_sil + ice.free.days, random = ~ stn 
                     + spl(par) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days) + stn:ma(profile.depth)
                     , na.method.X = " include ", na.method.Y = "include", 
                     data = glm.full, maxiter = 30, workspace = 400000000)


#------------------------------------- 12. RANDOM AR(2) DEPTH TERM ------------------------------#


fluoro.mm12 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                     + bot_sil + ice.free.days, random = ~ stn 
                     + spl(par) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days) + stn:ar2(profile.depth)
                     , na.method.X = " include ", na.method.Y = "include", 
                     data = glm.full.cut, maxiter = 30, workspace = 400000000)



#----------------------- 13. MIXED MODEL WITH RANDOM ARMA DEPTH TERM ----------------------------#


fluoro.mm13 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                      + bot_sil + ice.free.days, random = ~ stn 
                      + spl(par) + spl(temp) + spl(bot_oxy) 
                      + spl(bot_sil) + spl(ice.free.days) + stn:arma(profile.depth)
                      , na.method.X = " include ", na.method.Y = "include", 
                      data = glm.full.cut, maxiter = 30, workspace = 400000000)

#--------------------- 14. GAUSSIAN DEPTH ERROR STRUCTURE -------------------------#


fluoro.mm14 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                     + bot_sil + ice.free.days, random = ~ stn 
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                     + spl(bot_sil) + spl(ice.free.days) 
                     , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:gau(profile.depth),
                     data = glm.full, maxiter = 30, workspace = 400000000)


#---------------------------------- MIXED MODEL PLOTS -----------------------------------------#

fluoro.mm <- fluoro.mm1


summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)

#plot residuals for a certain station
s <- 40
plot(as.numeric(glm.full$profile.depth[glm.full$stn == s]), 
     fluoro.mm$residuals[glm.full$stn == s], xlab = "profile depth (m)", ylab = "residuals")
title(bquote("Plot of residuals at station" ~ .(s)~""))


#plot the residuals at 40m as a bubble plot by latitude and longitude
d <- 40
l.fluoro.bubble <- fluoro.mm$residuals[glm.full$profile.depth == d]
rad   <- sqrt(abs(l.fluoro.bubble)/pi) #set size = radius rather than area of circle
lat   <- glm.full$stn.lat.full[glm.full$profile.depth == d]
long  <- glm.full$stn.long.full[glm.full$profile.depth == d]
point_col <- rep(0, length(l.fluoro.bubble))
point_col[l.fluoro.bubble <= 0] <- "#0000FF"
point_col[l.fluoro.bubble > 0] <- "#FF0000"
symbols(x = long, y = lat, circles = rad, inches = 0.2, fg = point_col, lwd = 2.5)
title(bquote("Residuals by latitude/longitude at " ~ .(d)~"m depth"))
legend("topleft", c("negative", "positive"), col = c("blue", "red"), lwd = 2.5, bty = "n", cex = 0.75)


#correlogram to assess spatial auto-correlation
d <- 40
correlog <- spline.correlog(x = glm.full$stn.lat.full[glm.full$profile.depth == d], 
                            y = glm.full$stn.long.full[glm.full$profile.depth == d], 
                            z = fluoro.mm$residuals[glm.full$profile.depth == d], na.rm = T)
plot.spline.correlog(correlog)
title(bquote("Correlogram at " ~ .(d)~"m depth"))


#variogram at depth 40
d <- 40
v1 <- variog(coords = cbind(glm.full$stn.lat.full[glm.full$profile.depth == d], 
                            glm.full$stn.lat.full[glm.full$profile.depth == d]), 
             data = fluoro.mm$residuals[glm.full$profile.depth == d])
plot(v1, main = bquote("Semivariogram at " ~ .(d)~"m depth"), pch = 19)







