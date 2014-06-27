#runs linear mixed effects model on thinned data set using CTD data only in lme4

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
glm.full <- read.csv(file = "rstnCTD.csv", header= T) 
library(lmeSplines)
library(lme4)
attach(glm.full)

#remove NA values and reduce the number of stations used
glm.full[glm.full$fluoro == -9, c(9:17)] <- NA
glm.full <- glm.full[glm.full$stn %in% seq(2, 118, by = 10), ] #remove stations to cut down data set

#gamma mixed model using glmer
fluoro.mm <- glmer((fluoro + 1)~ par + profile.depth + temp + oxy + (1 | stn), family = Gamma, data = glm.full)
summary(fluoro.mm)
plot(fluoro.mm)
hist(resid(fluoro.mm))

s <- 2
plot(glm.full$profile.depth[glm.full$stn == s], glm.full$fluoro[glm.full$stn == s] + 1, xlab = "profile.depth", ylab = "", pch = 19)
points(glm.full$profile.depth[glm.full$stn == s], fitted(fluoro.mm)[glm.full$stn == s], col = "red", pch = 19)
title("Fitted (red) vs log fluoro (black) values for all depths of station 2")



#gamma mixed model using glmer with random splines
fluoro.mm <- glmer((fluoro + 1)~ par + profile.depth + temp + oxy + (bs(par) | stn) + (bs(profile.depth) | stn) + (bs(temp) | stn)+ (bs(oxy) | stn), family = Gamma, data = glm.full)
summary(fluoro.mm)
plot(fluoro.mm)
hist(resid(fluoro.mm))

s <- 2
plot(glm.full$profile.depth[glm.full$stn == s], glm.full$fluoro[glm.full$stn == s] + 1, xlab = "profile.depth", ylab = "", pch = 19)
points(glm.full$profile.depth[glm.full$stn == s], fitted(fluoro.mm)[glm.full$stn == s], col = "red", pch = 19)
title("Fitted (red) vs log fluoro (black) values for all depths of station 2")
