dat <- read.csv(file = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data/test data/ctd_example.csv", header= T)
names(dat) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "depth", "transmittance", "cond", "temp", "sal", "par", "oxygen", "fluoro", "x2", "ice", "wm")

test_asremlAIC <- function() {
  
  asreml.fit <- asreml(fixed = fluoro ~ depth, data = dat, na.method.X = "include", trace = F)
  
  answer <- asremlAIC(asreml.fit)
  
  checkEqualsNumeric(935, round(answer$AIC))
  checkEqualsNumeric(1, answer$K)
  checkEqualsNumeric(-466, round(answer$l))
  
}

test_calcRsquared <- function() {
  
  asreml.fit <- asreml(fixed = fluoro ~ depth, random =~ stn, data = dat, na.method.X = "include", trace = F)
  
  answer <- calcRsquared(asreml.fit, varStruct = FALSE, rand = "stn", data = dat)
  
  checkEqualsNumeric(0.33, round(answer$marginal, 2))
  checkEqualsNumeric(0.33, round(answer$conditional, 2))

}

test_deg2rad <- function() {
  
  checkEqualsNumeric(pi, deg2rad(180))
  
}

test_depthFluoroMax <- function() {
  
  suppressWarnings(dat$l.fluoro <- log(dat$fluoro))

  answer <- depthFluoroMax(dat)
  
  checkEqualsNumeric(25, answer[4])
  
}

test_distFromStn1 <- function() {
  
  lat  <- dat$lat[duplicated(dat$stn) == FALSE]
  long <- dat$long[duplicated(dat$stn) == FALSE]
  
  answer <- distFromStn1(lat, long)
  
  checkEqualsNumeric(0.89, round(answer$x[1], 2))
  checkEqualsNumeric(0, answer$y[1])
  
}

test_gcdHF <- function() {
  
  source("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/deg2rad.R")
  
  answer <- gcdHF(deg2rad(lat[1]), deg2rad(long[1]), deg2rad(lat[2]), deg2rad(long[2]))
  
  checkEqualsNumeric(92, round(answer))
  
}









