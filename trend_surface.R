#runs mixed model on thinned data set using CTD data only

setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
glm.full <- read.csv(file = "thinCTD.csv", header= T)
library(asreml)
attach(glm.full)

l.fluoro <- log(glm.full$fluoro)
l.fluoro[is.nan(l.fluoro) == T] <- NA
l.fluoro[l.fluoro == -Inf] <- NA
glm.full <- cbind(glm.full, l.fluoro)

#standard mixed model with square root transformation
fluoro.mm1 <- asreml(fixed = sqrt(fluoro) ~ par + profile.depth + temp + oxy, random = ~ stn 
                     + spl(par) + spl(profile.depth) + spl(temp) + spl(oxy)  
                     , na.method.X = "include", na.method.Y = "include", 
                     data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm1)
wald(fluoro.mm1)
plot(fluoro.mm1)
plot(profile.depth[stn == 100], fluoro.mm1$residuals[stn == 100])

#mixed model with square root transformation and autocovariate
autocov <- 0
for(i in 1:nrow(glm.full)){
  autocov[i] <- sum(sqrt(fluoro[i - 1]), sqrt(fluoro[i + 1]))/2
}
  
fluoro.mm1 <- asreml(fixed = sqrt(fluoro) ~ par + profile.depth + temp + oxy + autocov, 
                     random = ~ stn + spl(par) + spl(profile.depth) + spl(temp) + spl(oxy)  
                     , na.method.X = "include", na.method.Y = "include", 
                     data = glm.full, maxiter = 30, workspace = 400000000)
summary(fluoro.mm1)
wald(fluoro.mm1)
plot(fluoro.mm1)
plot(profile.depth[stn == 50], fluoro.mm1$residuals[stn == 50])




