#runs linear mixed effects model on thinned data set using CTD data only

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
#usually use thinCTD.csv, fullCTD.csv is the full data set using CTD data only and padding missing depths with NA
glm.full <- read.csv(file = "rstnCTD.csv", header= T) 
library(nlme)
library(lmeSplines)
attach(glm.full)

glm.full[glm.full$fluoro == -9, c(9:17)] <- NA
glm.full <- glm.full[glm.full$stn %in% seq(2, 118, by = 10), ] #remove stations to cut down data set


#non linear mixed effects model without correlation structure

fluoro.mm1 <- nlme(l.fluoro ~ exp(par) + profile.depth + profile.depth^2 + temp + oxy + stn, fixed = par + profile.depth + temp + oxy ~ 1, 
                  data = glm.full, na.action = na.omit, start = rep(1, 4), groups =~ stn)
summary(fluoro.mm1)
hist(fluoro.mm1$resid)
plot(fluoro.mm1)
plot(fluoro.mm1$resid[stn == 38])


#lme fit with single random spline
glm.spl <- na.omit(glm.full)
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ par, data = glm.spl)

fluoro.mm1 <- lme(l.fluoro ~ par + profile.depth + temp + oxy + stn, random=list(all=pdIdent(~Zt - 1)), data = glm.spl)
summary(fluoro.mm1)
hist(fluoro.mm1$resid)
plot(fluoro.mm1)
plot(fluoro.mm1$resid[stn == 38])


#lme fit with multiple random splines
# The fit for each station for this model looks pretty good
glm.spl <- na.omit(glm.full)
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ par, data = glm.spl)
glm.spl$Ztd  <- smspline(~ profile.depth, data = glm.spl)
glm.spl$Zto  <- smspline(~ oxy, data = glm.spl) #slow
glm.spl$Ztt  <- smspline(~ temp, data = glm.spl) #slow

fluoro.mm1 <- lme(l.fluoro ~ par + profile.depth + temp + oxy + stn, random=list(all=pdIdent(~Zt - 1), all=pdIdent(~Ztd - 1), all=pdIdent(~Zto - 1), all=pdIdent(~Ztt - 1), all =~ 1|stn), data = glm.spl)
summary(fluoro.mm1)
hist(resid(fluoro.mm1), main = "Histogram of residuals")
plot(fluoro.mm1, pch = 19)
plot(resid(fluoro.mm1)[glm.spl$stn == 2], pch = 19)

s <- 2
plot(glm.spl$profile.depth[glm.spl$stn == s], glm.spl$l.fluoro[glm.spl$stn == s], xlab = "profile.depth", ylab = "", pch = 19)
points(glm.spl$profile.depth[glm.spl$stn == s], fitted(fluoro.mm1)[glm.spl$stn == s], col = "red", pch = 19)
title("Fitted (red) vs log fluoro (black) values for all depths of station 2")


fluoro.mm1 <- lme(l.fluoro ~ par + oxy + profile.depth + temp, random= list(stn= pdBlocked(list(~ profile.depth, pdIdent(~Ztd - 1)))), data = glm.spl)

#centre and scale data
glm.centre <- na.omit(glm.full)
glm.centre$l.fluoro <- (glm.centre$l.fluoro - mean(glm.spl$l.fluoro))/sd(glm.spl$l.fluoro)
glm.centre$temp <- (glm.centre$temp - mean(glm.spl$temp))/sd(glm.spl$temp)
glm.centre$profile.depth <- na.omit(glm.full)$profile.depth
glm.centre$par <- (glm.centre$par - mean(glm.spl$par))/sd(glm.spl$par)
glm.centre$oxy <- (glm.centre$oxy - mean(glm.spl$oxy))/sd(glm.spl$oxy)

#model run with centred and scaled data
glm.spl <- na.omit(glm.centre)
glm.spl$all <- rep(1, nrow(glm.spl))
glm.spl$Zt  <- smspline(~ par, data = glm.spl)
glm.spl$Ztd  <- smspline(~ profile.depth, data = glm.spl)
glm.spl$Zto  <- smspline(~ oxy, data = glm.spl) #slow
glm.spl$Ztt  <- smspline(~ temp, data = glm.spl) #slow
fluoro.mm1 <- lme(l.fluoro ~ par + profile.depth + temp + oxy + stn, random=list(all=pdIdent(~Zt - 1), all=pdIdent(~Ztd - 1)), data = glm.spl,
                  correlation = corExp(form =~ profile.depth|stn))



#lme with random splines and linear correlation structure for grouped data
glm.spl <- na.omit(glm.full)
glm.spl$all <- rep(1, nrow(glm.spl)) #no grouping
glm.spl$Zt  <- smspline(~ par, data = glm.spl)
glm.spl$Ztd  <- smspline(~ profile.depth, data = glm.spl)
glm.spl$Zto  <- smspline(~ oxy, data = glm.spl) #slow
glm.spl$Ztt  <- smspline(~ temp, data = glm.spl) #slow

#use pdBlocked structure by station to indicate the station level blocking for each variable
fluoro.mm1 <- lme(l.fluoro ~ par + profile.depth - 1, random = list(stn= pdIdent(~Ztd - 1), stn= pdIdent(~Zto - 1), stn= pdIdent(~Ztt - 1),
                    stn= pdIdent(~Zt - 1), stn =~1), data = glm.spl)
summary(fluoro.mm1)
hist(resid(fluoro.mm1))
plot(fluoro.mm1)
plot(resid(fluoro.mm1)[glm.spl$stn == 10])


s <- 10
plot(glm.spl$profile.depth[glm.spl$stn == s], glm.spl$l.fluoro[glm.spl$stn == s], xlab = "profile.depth")
points(glm.spl$profile.depth[glm.spl$stn == s], fitted(fluoro.mm1)[glm.spl$stn == s], col = "red")

z <- na.omit(l.fluoro)
zlen <- length(l.fluoro)
levels <- seq(min(z),max(z),length.out = zlen)
col <- colorRampPalette(c("black","red"))(zlen)[rank(z)]
plot(temp, par, col = col, pch = 19)

plot(temp, oxy, col = col, pch = 19)


fluoro.mm <- glmer((fluoro + 1)~ par + profile.depth + temp + oxy + (1 | stn), family = Gamma, data = glm.full)
summary(fluoro.mm)
plot(fluoro.mm)
hist(resid(fluoro.mm))




s <- 22
plot(glm.full$profile.depth[glm.full$stn == s], glm.full$fluoro[glm.full$stn == s] + 1, xlab = "profile.depth", ylab = "", pch = 19)
points(glm.full$profile.depth[glm.full$stn == s], fitted(fluoro.mm)[glm.full$stn == s], col = "red", pch = 19)
title("Fitted (red) vs log fluoro (black) values for all depths of station 2")



plot(fitted(fluoro.mm)[stn == 62])




