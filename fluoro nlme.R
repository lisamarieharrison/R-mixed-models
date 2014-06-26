setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
dat <- read.csv(file = "thinCTD.csv", header = T)
attach(dat)
library(nlme)
library(lme4)
library(mgcv)

glm.dat <- dat[, c("fluoro", "sal", "temp", "par", "oxy", "stn", "profile.depth")]
glm.dat$stn <- as.factor(glm.dat$stn)
glm.dat$sal[glm.dat$sal == -9] <- NA
glm.dat$temp[glm.dat$temp == -9] <- NA
glm.dat$fluoro[glm.dat$fluoro <= 0] <- NA


#glm
fluoro.glm <- glm(fluoro ~ sal + par + temp + oxy + profile.depth, data = na.omit(glm.dat), family = Gamma(link = "inverse")
                  , start = c(1, 1, 1, 1, 1, 1))
summary(fluoro.glm)
plot(fluoro.glm)


#lme
fluoro.lme <- lme(fluoro ~ sal + par + temp + oxy + profile.depth, random =~ profile.depth|stn, 
                  data = na.omit(glm.dat))
summary(fluoro.lme)
plot(fluoro.lme)

#lme with log(fluoro)
glm.dat$l.fluoro <- log(glm.dat$fluoro)

fluoro.lme <- lme(l.fluoro ~ sal + par + temp + oxy + profile.depth, random =~ stn, 
                  data = na.omit(glm.dat))
summary(fluoro.lme)
par(mfrow = c(2, 2))
hist(resid(fluoro.lme))
qqnorm(resid(fluoro.lme))
plot(fitted(fluoro.lme), resid(fluoro.lme), xlab = "fitted", ylab = "residuals")
plot(resid(fluoro.lme), ylab = "residuals")

knit2html("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/fluoro glms.Rmd", "test.doc")

#lme with log(fluoro) and different slope for each station
glm.dat$l.fluoro <- log(glm.dat$fluoro)

fluoro.lme <- lme(l.fluoro ~ sal + par + temp + oxy + profile.depth, random =~ profile.depth|stn, 
                  data = na.omit(glm.dat))
summary(fluoro.lme)
par(mfrow = c(2, 2))
hist(resid(fluoro.lme))
qqnorm(resid(fluoro.lme))
plot(fitted(fluoro.lme), resid(fluoro.lme), xlab = "fitted", ylab = "residuals")
plot(resid(fluoro.lme), ylab = "residuals")

#lme with log(fluoro) and correlation structure
glm.dat$l.fluoro <- log(glm.dat$fluoro)

fluoro.lme <- lme(l.fluoro ~ sal + par + temp + oxy + profile.depth, random =~ 1|profile.depth, 
                  correlation = corAR1(form =~ 1|stn), data = na.omit(glm.dat), control = list(msMaxIter = 1000000))
summary(fluoro.lme)
par(mfrow = c(2, 2))
hist(resid(fluoro.lme))
qqnorm(resid(fluoro.lme))
plot(fitted(fluoro.lme), resid(fluoro.lme), xlab = "fitted", ylab = "residuals")
plot(resid(fluoro.lme), ylab = "residuals")

#marginal model

fluoro.mar <- gls(l.fluoro ~ sal + par + temp + oxy + profile.depth, method = "REML", correlation 
                  = corAR1(form =~ 1|stn), data = na.omit(glm.dat))
summary(fluoro.mar)
par(mfrow = c(2, 2))
hist(resid(fluoro.mar))
qqnorm(resid(fluoro.mar))
plot(fitted(fluoro.mar), resid(fluoro.mar), xlab = "fitted", ylab = "residuals")
plot(resid(fluoro.mar), ylab = "residuals")
knit2html("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/fluoro glms.Rmd", "test.doc")


#lmer with random station effect
fluoro.lmer <- lmer(l.fluoro ~ sal + par + temp + oxy + profile.depth + (1|stn), data = na.omit(glm.dat))
summary(fluoro.lmer)
par(mfrow = c(2, 2))
hist(resid(fluoro.lmer))
qqnorm(resid(fluoro.lmer))
plot(fitted(fluoro.lmer), resid(fluoro.lmer), xlab = "fitted", ylab = "residuals")
plot(resid(fluoro.lmer), ylab = "residuals")
knit2html("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/fluoro glms.Rmd", "test.doc")

#lmer with random station intercept and slope




