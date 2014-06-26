a <- which(dat$Cast.Number == 1 | dat$Cast.Number ==119 | dat$Cast.Number == 120)

dat_stn <- dat[-a, ]

a <- which(dat_stn$Pressure %in% seq(2, 250, by = 10))

dat_thin <- dat_stn[a, ]


names(dat) <- c("survey", "stn", "latitude", "longitude", "start.time", "end.time", "depth", 
                "transmittance.null", "conductivity.null", "temp", "salinity", "par", "oxygen", 
                "fluoro", "x2.dbar")

dat$PAR.uE.m..2.sec..uE.m..2.sec[dat$PAR.uE.m..2.sec..uE.m..2.sec == -9] <- NA
