#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL MA DEPTH CORRELATION ------------#

fluoro.mm1 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days)
                    , na.method.X = " include ", na.method.Y = "include", rov =~ stn:ma(profile.depth), 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)


#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL EXP DEPTH CORRELATION ------------#

fluoro.mm2 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days)
                    , na.method.X = " include ", na.method.Y = "include", rcov =~ stn: exp(profile.depth), 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)

#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RANDOM MA DEPTH TERM ------------#

fluoro.mm3 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days) + stn:ma(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)


#------------ MIXED MODEL WITH RANDOM AR(2) DEPTH TERM ------------#


fluoro.mm4 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days) + stn:ar2(profile.depth)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.cut, maxiter = 30, workspace = 400000000)




#------------------- MIXED MODEL RUN WITHOUT TRANSFORMATION TO LOG FLUORO ----------------------#


fluoro.mm5 <- asreml(fixed = fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days) 
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)








#--------------------- MIXED MODEL WITH SQUARE ROOT TRANSFORMATION -------------------------#


fluoro.mm6 <- asreml(fixed = sqrt(fluoro) ~ par + profile.depth + bot_nit + sal + temp + bot_oxy 
                    + bot_sil + ice.free.days + spl(bot_phos), random = ~ stn 
                    + spl(bot_nit) + spl(par) + spl(bot_oxy) + spl(bot_sil)
                    + spl(temp) + spl(ice.free.days) +spl(sal) + spl(bot_phos)
                    , na.method.X = " include ", na.method.Y = "include", 
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)



#--------------------- MIXED MODEL WITH GAUSSIAN DEPTH ERROR STRUCTURE -------------------------#


fluoro.mm7 <- asreml(fixed = l.fluoro ~ par + profile.depth + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(profile.depth) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days) 
                    , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:gau(profile.depth),
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)




#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL AR(1) DEPTH CORRELATION ------------#

fluoro.mm8 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days) 
                    , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:ar1(profile.depth),
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)



#------------ MIXED MODEL WITH ONLY SIGNIFICANT VARIABLES + RESIDUAL AR(2) DEPTH CORRELATION ------------#

fluoro.mm9 <- asreml(fixed = l.fluoro ~ par + temp + bot_oxy 
                    + bot_sil + ice.free.days, random = ~ stn 
                    + spl(par) + spl(temp) + spl(bot_oxy) 
                    + spl(bot_sil) + spl(ice.free.days) 
                    , na.method.X = " include ", na.method.Y = "include", rcov =~ stn:ar2(profile.depth),
                    data = glm.full.sort, maxiter = 30, workspace = 400000000)





