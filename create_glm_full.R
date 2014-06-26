setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.data  <- read.csv(file = "glm_data.csv")
extra.glm <- read.csv("extra_glm_data_july_06.csv")
stn.info  <- as.matrix(read.table(file = "station_information.txt", skip = 1, fill = T))


glm.data <- cbind(glm.data, extra.glm)
glm.data$water_mass <- as.factor(glm.data$water_mass)
glm.data$stn <- as.factor(glm.data$stn)
glm.data$l.fluoro <- log(glm.data$fluoro)
glm.data$l.fluoro[is.nan(glm.data$l.fluoro) == T] <- NA #change nan to NA
glm.data$l.fluoro[glm.data$l.fluoro == -Inf] <- NA #change nan to NA
attach(glm.data)

## -------------- FILL IN MISSING DEPTH DATA WITH NAS ---------------------- ##

#create matrix of the correct size with all station and depths included
glm.full <- matrix(NA, nrow = 14750, ncol = ncol(glm.data))
glm.full[, 8] <- rep(seq(2, 250, by = 2), 118)
glm.full[, 2] <- sort(rep(c(1:118), 125))

#fill in rows that we have information about
for (i in 1:nrow(glm.data)){
  j <- which(glm.full[, 2] == glm.data[i, 2] & glm.full[, 8] == glm.data[i, 8])
  glm.full[j, ] <- as.numeric(glm.data[i, ])
  if (i %% 1000 == 0) print(i)
  flush.console()
}

#add day column for each station
glm.full <- cbind(glm.full, 0, 0)
for (i in 1:118){
  w <- which(glm.full$stn == i)
  glm.full[w, 28] <- stn.info[i, 3]
  glm.full[w, 29] <- stn.info[i, 4]
}
glm.full[, 28]<- as.numeric(glm.full[, 28])

#add a column for days after start of survey
glm.full <- cbind(glm.full, 0)
start_day <- glm.full[1, 28]
for (i in 1:nrow(glm.full)){
  day <- glm.full[i, 28]
  month <- glm.full[i, 29]
  if(month == "Jan") glm.full[i, 30] <- day - start_day
  if(month == "Feb") glm.full[i, 30] <- day + 31 - start_day
}

glm.full <- as.data.frame(glm.full) #asreml requires a data frame
names(glm.full) <- names(glm.data)
names(glm.full)[28:30] <- c("day", "month", "time")
attach(glm.full)