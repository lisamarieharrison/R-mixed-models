setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.full <- read.table(file = "glm.full.txt", header= T)
stn.info <- read.table(file = "station_information.txt", skip = 1, fill = T)
attach(glm.full)
library(asreml)
library(calibrate)
library(ncf)
library(nlme)


# ------------------------- EXPLORATORY PLOTS ----------------------------- #

#plot par against depth
plot(profile.depth, par)

#plot par against log fluoro
plot(par, l.fluoro)

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

#plot l.fluoro against latitude/longitude for 40m depth
d <- 40
rbPal <- colorRampPalette(c('red','blue'))
point_col <- na.omit(l.fluoro[profile.depth == d]) - min(na.omit(l.fluoro[profile.depth == d]))
dat_cols <- rbPal(length(point_col))[as.numeric(cut(point_col,breaks = length(point_col)))]
point_labels <- time[profile.depth == d & latitude < -55]


plot(c(longitude[profile.depth == d & latitude < -55], c(80:89)), c(latitude[profile.depth == d & latitude < -55], rep(NA, 10)), 
     col = dat_cols, xlab = "longitude", ylab = "latitude", pch = 16)
title("Log fluoro at 40m depth at each station labelled with days since start of survey")
textxy(longitude[profile.depth == d & latitude < -55], latitude[profile.depth == d & latitude < -55],
       point_labels)
legend("right", legend = seq(-1.5, 2, 0.5), bty = "n",
       col = rbPal(12.5), pch = 16, cex = 0.6)

#plot l.fluoro against latitude/longitude for 40m depth - bubble plot
l.fluoro.bubble <- l.fluoro[longitude < 85 & profile.depth == 40]
rad   <- sqrt(abs(resid)/pi) #set size = radius rather than area of circle
lat   <- latitude[profile.depth == 40 & longitude < 85]
long  <- longitude[profile.depth == 40 & longitude < 85]
symbols(x = long, y = lat, circles = rad, inches = 0.3)
title("Log fluoro for depth = 40m")

#correlogram to assess spatial auto-correlation at depth = 40m
d <- 150
correlog <- spline.correlog(x = glm.full$latitude[glm.full$profile.depth == d], 
                            y = glm.full$longitude[glm.full$profile.depth == d], 
                            z = glm.full$l.fluoro[glm.full$profile.depth == d], na.rm = T)
plot.spline.correlog(correlog)
main_title <- bquote("Correlogram at " ~ .(d)~"m depth")
title(main_title)

#plot of bottled nitrogen against log fluoro at 40m depth
plot(bot_nit[latitude < 85 & profile.depth == 40], l.fluoro[latitude < 85 & profile.depth == 40],
     xlab = "bottled nitrogen", ylab = "l.fluoro")
title("log fluoro against bottled nitrogen at 40m depth")

#check for spatial autocorrelation of residuals
d <- 40
resid <- fluoro.mm$residuals[profile.depth == d & longitude < 85] + 2
correlog <- spline.correlog(x = latitude[profile.depth == d & longitude < 85], 
                            y = longitude[profile.depth == d & longitude < 85], 
                            z = resid, na.rm = T)
plot.spline.correlog(correlog)
main_title <- bquote("Correlogram of residuals at " ~ .(d)~"m depth")
title(main_title)

# ----------------------------- MIXED MODELS ------------------------------ #

#mixed model with depth, par, salinity, temperature, nitrogen and station
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + sal + temp 
                    + bot_nit, random = ~ stn + spl(profile.depth)
                    + spl(bot_nit) + spl(sal) + spl(temp) + spl(par), 
                    na.method.X = "omit", na.method.Y = "include", 
                    data = glm.full)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)

#mixed model with all variables
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + sal + temp + bot_oxy + 
                      bot_sil + bot_phos + bot_nit, random = ~ stn + spl(profile.depth)
                    + spl(bot_nit) + spl(par) + spl(sal) + spl(temp) + spl(bot_oxy)
                    + spl(bot_sil) + spl(bot_phos), na.method.X = "omit", 
                    na.method.Y = "omit", data = glm.full)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)

#mixed model with all variables and lat/long
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + sal + temp + bot_oxy 
                    + bot_sil + bot_phos + bot_nit, random = ~ stn + latitude*longitude 
                    + spl(profile.depth)
                    + spl(bot_nit) + spl(par) + spl(sal) + spl(temp) + spl(bot_oxy)
                    + spl(bot_sil) + spl(bot_phos), na.method.X = "omit", na.method.Y = "omit", 
                    data = glm.full)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm) #anova produces the same results

#drop salinity
#mixed model with all variables and lat/long
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + temp + bot_oxy + bot_sil + bot_phos
                    + bot_nit, random = ~ stn + latitude*longitude + spl(profile.depth)
                    + spl(bot_nit) + spl(par) + spl(temp) + spl(bot_oxy) + spl(bot_sil) 
                    + spl(bot_phos), na.method.X = "include", 
                    na.method.Y = "include", data = glm.data)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#residuals plotted by lat/long
resid <- fluoro.mm$residuals[profile.depth == 40 & latitude < -50]
rad   <- sqrt(abs(resid)/pi)
lat   <- latitude[profile.depth == 40 & latitude < -50]
long  <- longitude[profile.depth == 40 & latitude < -50]
#change colour of circle edges
edge.col <- resid
edge.col[resid < 0] <- "red"
edge.col[resid == 0] <- "yellow"
edge.col[resid > 0] <- "blue"

symbols(x = long, y = lat, circles = rad, inches = 0.3, fg = edge.col)
title("Residuals for depth = 40m. red = negative, blue = positive")

plot(profile.depth[latitude < -50], fluoro.mm$residuals[latitude < -50], 
     xlab = "depth (m)", ylab = "residuals")
title("Residuals by depth")

plot(profile.depth[latitude < -50], l.fluoro[latitude < -50], xlab = "depth (m)"
     , ylab = "log fluoro", main = "log fluoro against depth for all stations")

#add NA to stations where not all depths are represented

#mixed model with all variables and exp ar process for depth
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + temp + bot_oxy + 
                    + bot_nit, random = ~ stn + spl(profile.depth)+ spl(bot_nit) 
                    + spl(par) + spl(temp) + spl(bot_oxy), data = glm.full, 
                    rcov= ~ stn:exp(profile.depth), na.method.X = "include", 
                    na.method.Y = "include", maxiter = 30, workspace = 40000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#mixed model with all variables and ar process for depth and latitude*longitude
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + temp + bot_oxy + 
                      + bot_nit, random = ~ stn + latitude*longitude + spl(profile.depth)+ spl(bot_nit) 
                    + spl(par) + spl(temp) + spl(bot_oxy), data = glm.full, 
                    rcov= ~ stn:exp(profile.depth), na.method.X = "include", 
                    na.method.Y = "include", maxiter = 30, workspace = 40000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)


#mixed model with all variables and ar process for depth and latitude*longitude, 
#includes ice free days, station depth and water mass
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + bot_oxy + water_mass + bot_nit, random = ~ stn + latitude*longitude 
                    + spl(profile.depth)+ spl(bot_nit) + spl(par) + spl(temp) + spl(water_mass) + spl(ice.free.days)
                    + spl(stn_depth) + spl(bot_oxy), data = glm.full, rcov= ~ stn:exp(profile.depth), 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#variogram at depth 40
library(geoR)
d <- 40
v1 <- variog(coords = cbind(na.omit(latitude[profile.depth == d]), na.omit(longitude[profile.depth == d])), 
             data = na.omit(l.fluoro[profile.depth == d]))
plot(v1, main = "Variogram of all stations at 40m depth")


#mixed model with all variables and ar process for depth and latitude*longitude,
#includes ice free days, station depth and water mass drop nitrogen
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + bot_oxy + water_mass, random = ~ stn + latitude*longitude 
                    + spl(profile.depth) + spl(par) + spl(temp) + spl(water_mass)
                    + spl(stn_depth) + spl(bot_oxy), data = glm.full, rcov= ~ stn:exp(profile.depth), 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 40000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#mixed model with all variables and ar process for depth and latitude*longitude, 
#includes ice free days, station depth and water mass
#drop nitrogen and oxygen
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + water_mass, random = ~ stn + latitude*longitude 
                    + spl(profile.depth) + spl(par) + spl(temp) + spl(water_mass)
                    + spl(stn_depth), data = glm.full, rcov= ~ stn:exp(profile.depth), 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 40000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#mixed model with all variables and ar process for depth and latitude*longitude, 
#includes ice free days, station depth and water mass
#drop nitrogen, oxygen and water mass
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth, 
                    random = ~ stn + latitude*longitude 
                    + spl(profile.depth) + spl(par) + spl(temp)
                    + spl(stn_depth), data = glm.full, rcov= ~ stn:exp(profile.depth), 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 40000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#mixed model with all variables and ar process for depth and latitude*longitude, 
#includes ice free days, station depth, water mass and station depth
#drop nitrogen, oxygen, water mass and station depth
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp, 
                    random = ~ stn + latitude*longitude 
                    + spl(profile.depth) + spl(par) + spl(temp), data = glm.full, 
                    rcov= ~ stn:exp(profile.depth), 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 40000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#mixed model with all variables and ar process for depth and latitude*longitude, 
#includes ice free days, station depth, water mass and station depth
#drop nitrogen, oxygen, water mass, station depth and ice free days
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + par + temp, 
                    random = ~ stn + latitude*longitude 
                    + spl(profile.depth) + spl(par) + spl(temp), data = glm.full, 
                    rcov= ~ stn:exp(profile.depth), 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 40000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#-------------------------------- INCLUDING SPATIAL CORRELATION -------------------------------------#


#mixed model with all variables and ar process for depth and latitude*longitude, 
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + bot_oxy + water_mass + bot_nit, random = ~ stn + latitude*longitude 
                    + spl(profile.depth)+ spl(bot_nit) + spl(par) + spl(temp) + spl(water_mass) + spl(ice.free.days)
                    + spl(stn_depth) + spl(bot_oxy), data = glm.full, rcov= ~ stn:exp(profile.depth), 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)


#removing latitude and longitude
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + bot_oxy + water_mass + bot_nit, random = ~ stn
                    + spl(profile.depth)+ spl(bot_nit) + spl(par) + spl(temp) + spl(water_mass) + spl(ice.free.days)
                    + spl(stn_depth) + spl(bot_oxy), data = glm.full, rcov= ~ stn:exp(profile.depth), 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#removing the depth autocorrelation term
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + bot_oxy + water_mass + bot_nit, random = ~ stn
                    + spl(profile.depth)+ spl(bot_nit) + spl(par) + spl(temp) + spl(water_mass) + spl(ice.free.days)
                    + spl(stn_depth) + spl(bot_oxy), data = glm.full, 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#removing the exponential component of the depth correlation term
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + bot_oxy + water_mass + bot_nit, random = ~ stn
                    + spl(profile.depth)+ spl(bot_nit) + spl(par) + spl(temp) + spl(water_mass) + spl(ice.free.days)
                    + spl(stn_depth) + spl(bot_oxy), rcov= ~ stn:profile.depth, data = glm.full, 
                    na.method.X = "include", na.method.Y = "include", 
                    maxiter = 30, workspace = 400000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#removing rcov and omitting nas'
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + bot_oxy + water_mass + bot_nit, random = ~ stn + latitude*longitude
                    + spl(profile.depth)+ spl(bot_nit) + spl(par) + spl(temp) + spl(water_mass) + spl(ice.free.days)
                    + spl(stn_depth) + spl(bot_oxy), data = glm.full, 
                    na.method.X = "omit", na.method.Y = "omit", 
                    maxiter = 30, workspace = 400000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)

#removing rcov and omitting nas' without latitude and longitude
fluoro.mm <- asreml(fixed = l.fluoro ~ profile.depth + ice.free.days + par + temp + stn_depth
                    + bot_oxy + water_mass + bot_nit, random = ~ stn
                    + spl(profile.depth)+ spl(bot_nit) + spl(par) + spl(temp) + spl(water_mass) + spl(ice.free.days)
                    + spl(stn_depth) + spl(bot_oxy), data = glm.full, 
                    na.method.X = "omit", na.method.Y = "omit", 
                    maxiter = 30, workspace = 400000000)
fluoro.mm <- update(fluoro.mm) #run a few more iterations so LL converges
summary(fluoro.mm)
plot(fluoro.mm)
wald(fluoro.mm)





