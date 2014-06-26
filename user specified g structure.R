
depth.cor <- corSymm(l.fluoro, profile.depth|stn)


#g structure
fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + bot_oxy + bot_sil + ice.free.days
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) + stn:depth.cor
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



#-------------------------------------- R STRUCTURE -------------------------------------------#

depth.cor <- corSymm(fluoro.mm$resid, ~ profile.depth|stn)
depth.cor <- corSymm(~ 1 | glm.na.omit$profile.depth/glm.na.omit$stn)

glm.na.omit <- na.omit(glm.full.sort)

b <- Initialize(depth.cor, glm.na.omit$l.fluoro)

glm.full.cor <- cbind(glm.full.sort, a)
glm.cor.sort <- glm.full.cor[order(a), ]

fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + bot_oxy + bot_sil + ice.free.days + sal
                    , random = ~ stn + stn.lat.full 
                    + stn.long.full + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil) + spl(sal)
                    + spl(temp) + spl(ice.free.days)) + 
                    , na.method.X = " omit ", na.method.Y = "omit", 
                    data = glm.cor.sort, maxiter = 30, workspace = 400000000)
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





fluoro.mm <- asreml(fixed = l.fluoro ~ par + bot_nit + temp + bot_oxy + bot_sil + ice.free.days + sal
                    , random = ~ stn 
                    + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil) + spl(sal)
                    + spl(temp) + spl(ice.free.days) 
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)
summary(fluoro.mm)
wald(fluoro.mm)
plot(fluoro.mm)




