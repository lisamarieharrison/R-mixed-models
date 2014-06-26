

#plot of station number by latitude/longitude 
plot(stn.long, stn.lat, col = "white", xlab = "longitude",
     ylab = "latitude", main = "station number by latitude/longitude")
text(stn.long, stn.lat, labels = c(1:118))
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
v1 <- variog(coords = cbind(na.omit(glm.full.sort)$stn.lat.full[na.omit(glm.full.sort)$profile.depth == d], 
                            na.omit(glm.full.sort)$stn.lat.full[na.omit(glm.full.sort)$profile.depth == d]), 
             data = na.omit(glm.full.sort)$l.fluoro[na.omit(glm.full.sort)$profile.depth == d])
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

