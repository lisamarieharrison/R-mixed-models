setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/data")
glm.full <- read.csv(file = "glm_full.csv", header= T)
attach(glm.full)


stn_coords <- cbind(glm.full$stn, glm.full$stn.lat.full, glm.full$stn.long.full)
stn_coords <- stn_coords[seq(1, 14625, by = 125), ]
names(stn_coords) <- c("station", "latitude", "longitude")

write.csv(stn_coords, file = "stn_coordinates.csv", row.names = F)
