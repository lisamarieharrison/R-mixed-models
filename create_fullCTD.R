#Creates the data set fullCTD.csv, which contains all stations and depths using only CTD data
#missing depths are padded with NA so that correlation structures can be fitted in asreml


setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
glm.full <- read.csv(file = "procCTD.csv", header= T) #full data set with missing depths
glm.full <- glm.full[glm.full[, 2] %in% c(2:118), ] #remove test stations


#pad data set with NA if depth is missing
mat <- matrix(NA, nrow = 14625, ncol = 15)
mat[, 2] <- sort(rep(c(2:118), 125))
mat[, 7] <- rep(seq(2, 250, 2), 117)
for(i in 1:nrow(mat)){
  s <- mat[i, 2]
  d <- mat[i, 7]
  mat[i, ] <- as.numeric(glm.full[which(glm.full[, 2] == s & glm.full[, 7] == d), ])
}
mat[, 2] <- sort(rep(c(2:118), 125))
mat[, 7] <- rep(seq(2, 250, 2), 117)
glm.full <- data.frame(mat)

#add column names
names(glm.full) <- c("survey", "stn", "lat", "long", "start.time", "end.time", "profile.depth", "trans", "cond.null", "temp", "sall", "par", "oxy", "fluoro", "x2.dbar")
attach(glm.full)

#calculate log fluoro and add to data frame
l.fluoro <- log(glm.full$fluoro)
l.fluoro[is.nan(l.fluoro) == T] <- NA
l.fluoro[l.fluoro == -Inf] <- NA
glm.full <- cbind(glm.full, l.fluoro)

#write.csv(glm.full, "fullCTD.csv")