#runs additive mixed models in R using mgcv'
setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
glm.full <- read.csv(file = "rstnCTD.csv", header= T) 
library(lmeSplines)
library(mgcv)
attach(glm.full)

#remove NA values and reduce the number of stations used
glm.full[glm.full$fluoro == -9, c(9:17)] <- NA
glm.full <- glm.full[glm.full$stn %in% seq(2, 118, by = 10), ] #remove stations to cut down data set


#additive model using gamm (no random effects)
fluoro.mm <- gamm((1 + fluoro )~ s(par) + s(profile.depth) + s(temp) + s(oxy), family = Gamma, data = glm.full)
summary(fluoro.mm$lme)
plot(fluoro.mm$gam)
gam.check(fluoro.mm$gam)
res <- residuals(fluoro.mm$gam) #extract residuals
fit <- fitted(fluoro.mm$gam) #extract residuals

s <- 2
plot(glm.full$profile.depth[glm.full$stn == s], glm.full$fluoro[glm.full$stn == s] + 1, xlab = "profile.depth", ylab = "", pch = 19)
points(glm.full$profile.depth[glm.full$stn == s], fit[glm.full$stn == s], col = "red", pch = 19)
title("Fitted (red) vs log fluoro (black) values for all depths of station 2")



#additive mixed model using gamm
fluoro.mm <- gamm((1 + fluoro )~ s(par) + s(profile.depth) + s(temp) + s(oxy), random = list(stn =~ 1), family = Gamma, data = glm.full)
summary(fluoro.mm$lme)
plot(fluoro.mm$gam)
gam.check(fluoro.mm$gam)
res <- residuals(fluoro.mm$gam) #extract residuals
fit <- fitted(fluoro.mm$gam) #extract residuals

s <- 2
plot(glm.full$profile.depth[glm.full$stn == s], glm.full$fluoro[glm.full$stn == s] + 1, xlab = "profile.depth", ylab = "", pch = 19)
points(glm.full$profile.depth[glm.full$stn == s], fit[glm.full$stn == s], col = "red", pch = 19)
title("Fitted (red) vs log fluoro (black) values for all depths of station 2")


#additive mixed model with correlated errors 
fluoro.mm <- gamm((1 + fluoro )~ s(par) + s(profile.depth) + s(temp) + s(oxy), random = list(stn =~ 1), 
                  family = Gamma, data = glm.full, correlation = corExp())
summary(fluoro.mm$lme)
plot(fluoro.mm$gam)
gam.check(fluoro.mm$gam)
res <- residuals(fluoro.mm$gam) #extract residuals
fit <- fitted(fluoro.mm$gam) #extract residuals

s <- 2
plot(glm.full$profile.depth[glm.full$stn == s], glm.full$fluoro[glm.full$stn == s] + 1, xlab = "profile.depth", ylab = "", pch = 19)
points(glm.full$profile.depth[glm.full$stn == s], fit[glm.full$stn == s], col = "red", pch = 19)
title("Fitted (red) vs log fluoro (black) values for all depths of station 2")

