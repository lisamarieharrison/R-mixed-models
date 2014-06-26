setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.full <- read.csv(file = "procCTD.csv", header= T)

glm.stn2 <- subset(glm.full, Cast.Number == 2)
names(glm.stn2) <- c("survey", "stn", "latitude", "longitude", "start.time", "end.time", "depth", 
                "transmittance.null", "conductivity.null", "temp", "salinity", "par", "oxygen", 
                "fluoro", "x2.dbar")
glm.stn2$fluoro[glm.stn2$fluoro == -9] <- NA

fluoro <- matrix(glm.stn2$fluoro, nrow = 1)
write.table(fluoro, "C:/Users/Lisa/Documents/phd/southern ocean/WinBUGS/data/fluoro.txt", row.names = F, col.names = F, sep  = ",")


depth <- matrix(glm.stn2$depth, nrow = 1)
write.table(depth, "C:/Users/Lisa/Documents/phd/southern ocean/WinBUGS/data/depth.txt", row.names = F, col.names = F, sep  = ",")

temp <- matrix(glm.stn2$temp, nrow = 1)
write.table(temp, "C:/Users/Lisa/Documents/phd/southern ocean/WinBUGS/data/temp.txt", row.names = F, col.names = F, sep  = ",")




