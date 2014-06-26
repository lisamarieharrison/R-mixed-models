setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models")
glm.full <- read.csv(file = "glm_full.csv", header= T)
attach(glm.full)


plot(profile.depth, glm.full$fluoro, main = "Fluoro vs profile depth for all stations", type = "l")
plot(profile.depth, l.fluoro, main = "Fluoro vs profile depth for all stations", type = "l")


plot(profile.depth, par, main = "Depth vs par")
plot(profile.depth, temp, main = "Depth vs temperature", col = water_mass)
plot(profile.depth, sal, main = "Depth vs salinity")
plot(profile.depth, bot_nit, main = "Depth vs bottled nitrogen")
plot(profile.depth, bot_sil, main = "Depth vs bottled silicon")
plot(profile.depth, bot_oxy, main = "Depth vs bottled oxygen")

plot(bot_nit, bot_sil, main = "Nitrogen vs silicon")
plot(bot_nit, bot_oxy, main = "Nitrogen vs oxygen")

#fluoro vs depth
plot(profile.depth, glm.full$fluoro, main = "Fluoro vs profile depth for all stations", pch = 19,
     col = water_mass)

#Truncated to remove higher values
plot(profile.depth[glm.full$fluoro < 8], glm.full$fluoro[glm.full$fluoro < 8], main = "Fluoro vs profile depth for all stations", 
     col = water_mass, pch = 19)

stn[fluoro > 8]

plot(profile.depth, temp)
points(profile.depth[stn == 96], temp[stn == 96], col = "red")



raw_stn3<- read.csv(file = "C:/Users/Lisa/Documents/phd/southern ocean/BROKE-West_CTD_data/Processed/CTDData_200506031_080704/200506030-3.csv", header = T)

plot(profile.depth[stn == 3], bot_oxy[stn == 3])
points(profile.depth[stn ==3], raw_stn3$Oxygen..micromole.litre..micromole.litre[raw_stn3$Pressure <= 250], col = "red")


plot(bot_oxy, bot_sil)
