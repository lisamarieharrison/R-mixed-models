# Analysis for Guy Williams (see glm_readme.doc)

library(lattice)
library(asreml)

# model 'non-zero' values of fluoro 

rm(glm.data)

glm.data <- read.csv(file="C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/glm_data.csv")

summary(glm.data)

names(glm.data)

extra.glm.data.july.06 <- read.csv("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/extra_glm_data_july_06.csv")

names(extra.glm.data.july.06)
dim(extra.glm.data.july.06)
dim(glm.data)

glm.data$ice.fdays <- as.double(extra.glm.data.july.06$ice.free.days)

glm.data$HMLZ <-as.factor(extra.glm.data.july.06$horizontal.mixed.layer.zone)

glm.data$HMLZ <- as.integer(glm.data$HMLZ)

glm.data$profile.depth.f <- as.factor(glm.data$profile.depth)

glm.data$stn.f <- as.factor(glm.data$stn)
names(glm.data)

Nstat <- length(levels(glm.data$stn.f))
Ndepth <- length(levels(glm.data$profile.depth.f))

Nobs <- dim(glm.data)[1]

Ndepth*Nstat

summary(glm.data$HMLZ)

Nobs.stn <- tapply(X=rep(1,Nobs), INDEX=glm.data$stn.f, FUN=sum)

stn.mv <- levels(glm.data$stn.f)[Nobs.stn<125]
data.temp <- glm.data[glm.data$stn.f==stn.mv[1],]
dim(data.temp)

ind.p <- outer(data.temp$profile.depth,as.integer(levels(glm.data$profile.depth.f)),FUN="==")%*%
      matrix(data=seq(1,125,1), nrow=125, ncol=1)

Nlen <- length(ind.p)

N.extra <- 125-Nlen

data.temp.m <- matrix(data=rep(data.temp[ind.p[length(ind.p)],],each=N.extra), nrow=N.extra, ncol=dim(data.temp)[2])
data.temp <- as.data.frame(data.temp.m[,1:26])
glm.data.e <- glm.data[,1:26]
names(data.temp) <- names(glm.data.e)
data.temp$fluoro <- NA
data.temp$profile.depth <- as.integer(levels(glm.data$profile.depth.f)[(Nlen+1):125])

glm.data.e <- rbind(glm.data.e,data.temp)

data.temp <- glm.data[glm.data$stn.f==stn.mv[2],]
dim(data.temp)

ind.p <- outer(data.temp$profile.depth,as.integer(levels(glm.data$profile.depth.f)),FUN="==")%*%
      matrix(data=seq(1,125,1), nrow=125, ncol=1)

Nlen <- length(ind.p)

N.extra <- 125-Nlen
data.temp.m <- matrix(data=rep(data.temp[ind.p[length(ind.p)],],each=N.extra), nrow=N.extra, ncol=dim(data.temp)[2])
data.temp <- as.data.frame(data.temp.m[,1:26])
names(data.temp) <- names(glm.data.e)
data.temp$fluoro <- NA
data.temp$profile.depth <- as.integer(levels(glm.data$profile.depth.f)[(Nlen+1):125])

glm.data.e <- rbind(glm.data.e,data.temp)

data.temp <- glm.data[glm.data$stn.f==stn.mv[3],]
dim(data.temp)

ind.p <- outer(data.temp$profile.depth,as.integer(levels(glm.data$profile.depth.f)),FUN="==")%*%
      matrix(data=seq(1,125,1), nrow=125, ncol=1)

Nlen <- length(ind.p)

N.extra <- 125-Nlen
data.temp.m <- matrix(data=rep(data.temp[ind.p[length(ind.p)],],each=N.extra), nrow=N.extra, ncol=dim(data.temp)[2])
data.temp <- as.data.frame(data.temp.m[,1:26])
names(data.temp) <- names(glm.data.e)
data.temp$fluoro <- NA
data.temp$profile.depth <- as.integer(levels(glm.data$profile.depth.f)[(Nlen+1):125])

glm.data.e <- rbind(glm.data.e,data.temp)

dim(glm.data.e)

#summary(glm.data.e$HMLZ)
glm.data.e$HMLZ <- as.integer(glm.data.e$HMLZ)

names(glm.data.e)
names(glm.data)
rm(glm.data)

glm.data.e$stn <- as.integer(glm.data.e$stn)
summary(glm.data.e$stn)
glm.data <- glm.data.e[order(glm.data.e$stn,glm.data.e$profile.depth),]


glm.data$HMLZ.f <-as.factor(glm.data$HMLZ)

glm.data$profile.depth.f <- as.factor(glm.data$profile.depth)

glm.data$stn.f <- as.factor(glm.data$stn)

levels(glm.data$HMLZ.f)
names(glm.data)
glm.data <- glm.data[,c((1:25),27:29)]

levels(glm.data$stn.f)
length(levels(glm.data$stn.f))

#summary(glm.data)
names(glm.data)



# fit LMM using asreml()



tapply(X=as.integer(glm.data$transect), INDEX=list(glm.data$stn.f), FUN=mean)

glm.data$stn.wt.f <- as.factor(as.double(glm.data$transect)*1000+as.double(glm.data$stn))

length(levels(glm.data$stn.wt.f))
length(levels(glm.data$stn.f))


# correct HMLZ.f 

glm.data$latitude <- unlist(glm.data$latitude)
glm.data$HMLZ.fc <- as.factor(as.integer(glm.data$HMLZ.f)+ (as.integer(glm.data$HMLZ.f)==1 & 
    (abs(glm.data$latitude)<64)))

tapply(X=rep(1,dim(glm.data)[1]), INDEX=glm.data$HMLZ.f, FUN=sum)
sum(as.integer(glm.data$HMLZ.f))
summary(abs(glm.data$latitude))

glm.data$L.fluoro <- log(glm.data$fluoro+1)

summary(glm.data$L.fluoro)
summary(glm.data$fluoro)

# plot fluoro vs depth for stations 1:19

glm.data.prof <- glm.data[glm.data$stn>85 & glm.data$stn < 103,]
glm.data.prof$stn.f <- as.factor(glm.data.prof$stn)
names(glm.data.prof)

glm.data.prof <- glm.data[glm.data$stn>25 & glm.data$stn < 45,]
glm.data.prof$stn.f <- as.factor(glm.data.prof$stn)
names(glm.data.prof)


tapply(X=glm.data$fluoro[!is.na(glm.data$fluoro)], INDEX=glm.data$stn.f[!is.na(glm.data$fluoro)], FUN=max)

max(as.vector(tapply(X=glm.data$fluoro[!is.na(glm.data$fluoro)], INDEX=glm.data$stn.f[!is.na(glm.data$fluoro)], FUN=max)))
glm.data$longitude <- unlist(glm.data$longitude)
lat.stn <- as.vector(tapply(X=glm.data$latitude, INDEX=glm.data$stn.f, FUN=max))
long.stn <- as.vector(tapply(X=glm.data$longitude, INDEX=glm.data$stn.f, FUN=max))

glm.data.prof$longitude <- unlist(glm.data.prof$longitude)
lat.stn <- as.vector(tapply(X=glm.data.prof$latitude, INDEX=glm.data.prof$stn.f, FUN=max))
long.stn <- as.vector(tapply(X=glm.data.prof$longitude, INDEX=glm.data.prof$stn.f, FUN=max))



summary(glm.data$longitude)

cbind(lat.stn,long.stn,as.integer(levels(glm.data$stn.f)))

plot(y=lat.stn, x=long.stn, pch=levels(glm.data.prof$stn.f),ylim=c(-70,-60), xlim=c(29,81), cex=0.5)

xyplot(profile.depth ~ fluoro  | stn.f, type="l", data=glm.data.prof,
      ylim=c(250,-0.5), scales=list(y=list(at=seq(250,0,-50), labels=c("250","200","150","100","50","0"))),
        layout=c(4,5,1))

# group ice.fdays

glm.data$ice.fdays <- unlist(glm.data$ice.fdays)

hist(glm.data$ice.fdays)

# fit log-linear model
names(glm.data)

glm.data$water.mass <- unlist(glm.data$water_mass)
glm.data$water.mass[1:10]

glm.data$water.mass.f <- as.factor(glm.data$water.mass)

levels(glm.data$water.mass.f)
glm.data$par <- unlist(glm.data$par)
glm.data$par[1:10]
summary(glm.data$ice.fdays)
ice.fdays.f <- factor(cut(glm.data$ice.fdays, breaks=seq(-40,110,10)))
ice.fdays.m <- tapply(X=glm.data$ice.fdays, INDEX=list(ice.fdays.f), FUN=mean)
levels(ice.fdays.f)
ice.fdays.ml <- rep(0,length(glm.data$ice.fdays))
ice.fdays.ml <- ice.fdays.m [ice.fdays.f]
ice.fdays.ml <- as.double(ice.fdays.ml)
glm.data$ice.fdays.ml <- ice.fdays.ml
rm(ice.fdays.f,ice.fdays.ml)
glm.data$bot_phos <- unlist(glm.data$bot_phos)
summary(glm.data$bot_phos)
summary(glm.data$ice.fdays.ml)

glm.data$temp <- unlist(glm.data$temp)
hist(glm.data$par[glm.data$par>0])
plot(y=glm.data$par, x=glm.data$profile.depth)
summary(glm.data$par)


par.f <- factor(cut(glm.data$par, breaks=c(seq(0.0,0.5,0.1),seq(0.6,2.4,0.4))))
par.m <- tapply(X=glm.data$par, INDEX=list(par.f), FUN=mean)
par.ml <- rep(0,length(glm.data$par))
glm.data$par.ml <- par.m [par.f]
temp.f <- factor(cut(glm.data$temp, breaks=seq(-2,11,0.5)))
temp.m <- tapply(X=glm.data$temp, INDEX=list(temp.f), FUN=mean)
glm.data$temp.ml <- rep(0,length(glm.data$temp))
glm.data$temp.ml <- temp.m [temp.f]

summary(glm.data$par.ml)

names(glm.data)

asreml.02<-asreml(fixed = L.fluoro ~ water.mass.f+profile.depth+par.ml,
  random = ~ spl(profile.depth)+spl(par.ml)+stn.f ,  data=glm.data,
  na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)

asreml.02<-asreml(fixed = L.fluoro ~ water.mass.f+profile.depth+par.ml,
  random = ~ spl(profile.depth)+spl(par.ml)+stn.f, rcov= ~ stn.f:exp(profile.depth),data=glm.data,
  na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)


asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+water.mass.f:ice.fdays.ml+profile.depth+par.ml,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ stn.f,
  rcov= ~ stn.f:exp(profile.depth),data=glm.data,
  na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)


asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+water.mass.f:ice.fdays.ml+profile.depth+par.ml.ml+temp.ml+bot_phos,
  random = ~ spl(profile.depth)+spl(par.ml.ml)+ spl(ice.fdays.ml)+ water.mass.f:spl(ice.fdays.ml) + stn.f,
  rcov= ~ stn.f:exp(profile.depth),data=glm.data,
  na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)

asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+water.mass.f:ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ water.mass.f:spl(ice.fdays.ml) + stn.f,
  data=glm.data, na.method.Y = "include", na.method.X = "include", subset=fluoro > 0, maxiter=30, workspace=4000000)

asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+water.mass.f:ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ water.mass.f:spl(ice.fdays.ml) + stn.f,
  data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)


asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+water.mass.f:ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ water.mass.f:spl(ice.fdays.ml) + stn.f,
  rcov= ~ stn.f:exp(profile.depth), data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)


glm.data$water.mass.ff <- as.integer(glm.data$water.mass.f)*(as.integer(glm.data$water.mass.f)<3)+3*(as.integer(glm.data$water.mass.f)>2)
glm.data$water.mass.ff <- as.factor(glm.data$water.mass.ff)
levels(glm.data$water.mass.ff)
names(glm.data)
levels(glm.data$HMLZ.fc)

aasreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ stn.f:diag(water.mass.ff),
  rcov= ~ stn.f:exp(profile.depth), data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)

asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos + water.mass.ff:HMLZ.fc,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ stn.f:us(water.mass.ff),
  rcov= ~ stn.f:exp(profile.depth), data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)


asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos + HMLZ.fc + water.mass.ff:HMLZ.fc,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ stn.f:diag(water.mass.ff),
  rcov= ~ stn.f:exp(profile.depth), data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)


asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos + HMLZ.fc,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ stn.f:diag(water.mass.ff),
  rcov= ~ stn.f:exp(profile.depth), data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)


asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ stn.f:diag(water.mass.f),
  rcov= ~ stn.f:exp(profile.depth), data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)


asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+profile.depth+par.ml+temp.ml+bot_phos,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(ice.fdays.ml)+ stn.f:diag(water.mass.ff),
  rcov= ~ stn.f:exp(profile.depth), data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)

asreml.02<-asreml(fixed = L.fluoro ~ ice.fdays.ml+profile.depth+par.ml+temp.ml,
  random = ~ spl(profile.depth)+spl(par.ml)+ spl(temp.ml) + spl(ice.fdays.ml)+ stn.f:diag(water.mass.f),
  rcov= ~ stn.f:exp(profile.depth), data=glm.data, na.method.Y = "include", na.method.X = "include", maxiter=30, workspace=4000000)



levels(glm.data$HMLZ.fc)
asreml.02$loglik

anova(asreml.02)

(summary(asreml.02))$varcomp

var.v <- ((summary(asreml.02))$varcomp)[,2]

#plot variograms by water.mass

lag <- rep(0:50, times=3)
phi <- var.v[9]
#phi <- 0.8
wm.fac <- factor(x=rep(1:3,each=51), levels=c(1:3), labels=c("wm1","wm2","wm3&4"))
vario <- var.v[5]*as.integer(wm.fac=="wm1") + var.v[6]*as.integer(wm.fac=="wm2") + var.v[7]*as.integer(wm.fac=="wm3&4") +
       var.v[8]*(1-phi^lag)
data.vario <- data.frame(lag,vario,wm.fac)

xyplot(vario ~ lag | wm.fac, type="l", data=data.vario, xlab="Lag", ylab="semi-variance", layout=c(3,1))


(summary(asreml.02))$coef.fixed

re <- (summary(asreml.02))$coef.random

#edit(re)

re.stn <- as.vector(re[69:422,1])

re.stn.wm <- matrix(data=re.stn, nrow = length(re.stn)/3, ncol=3, byrow=T)

re.stn[1:6]
re.stn.wm[1:2,]

hist(re.stn.wm[,1])
hist(re.stn.wm[,2])
hist(re.stn.wm[,3])

re.stn.wm1 <- re.stn.wm[,1]
re.stn.wm2 <- re.stn.wm[,2]
re.stn.wm3 <- re.stn.wm[,3]
stn.id <- as.integer(levels(glm.data$stn.f))

lat.stn <- as.vector(tapply(X=glm.data$latitude, INDEX=glm.data$stn.f, FUN=max))
long.stn <- as.vector(tapply(X=glm.data$longitude, INDEX=glm.data$stn.f, FUN=max))
HMLZ.stn <- as.vector(tapply(X=as.integer(glm.data$HMLZ.fc), INDEX=glm.data$stn.f, FUN=mean))

plot(y=lat.stn, x=long.stn, pch=HMLZ.stn, ylim=c(-70,-60), xlim=c(28,82))

plot(y=lat.stn, x=long.stn, pch=2, ylim=c(-70,-60), xlim=c(28,82))



lat.stn,long.stn
plot(y=re.stn.wm[,1], x=re.stn.wm[,3])

re.stn.wm[1:20,]


res.asreml.02 <- asreml.02$residuals

length(res.asreml.02)

fv.asreml.02 <- exp(asreml.02$fitted.values)-1

# calculated marginal preds

re.byu <- re.stn.wm1[glm.data$stn.f]*as.integer(glm.data$water.mass.ff=="1")+
        re.stn.wm2[glm.data$stn.f]*as.integer(glm.data$water.mass.ff=="2") +
        re.stn.wm3[glm.data$stn.f]*as.integer(glm.data$water.mass.ff=="3")

summary(re.byu)

fv.asreml.02 <- exp(asreml.02$fitted.values + re.byu)-1
var(fv.asreml.02[!is.na(fv.asreml.02)])^0.5

summary(fv.asreml.02)

summary(glm.data$L.fluoro)

fv.asreml.02 <- exp(asreml.02$fitted.values)-1

fv.asreml.02 <- exp(asreml.02$fitted.values + re.byu)-1


index.sel <- glm.data$water.mass.f=="2" & glm.data$par < 0.1

summary(fv.asreml.02[index.sel])


fv.asreml.02.stn <- as.vector(tapply(X=fv.asreml.02[index.sel], INDEX=glm.data$stn.f[index.sel], FUN=max))

fluoro.stn <- as.vector(tapply(X=glm.data$fluoro[index.sel], INDEX=glm.data$stn.f[index.sel], FUN=max))


plot(x=fv.asreml.02[fv.asreml.02>0], y=glm.data$fluoro[fv.asreml.02>0], xlab="Conditional Prediction", ylab="Relative Fluoro",
       cex=0.4)
Pfluoro <- fv.asreml.02[fv.asreml.02>0]
Pfluoro <- Pfluoro[!is.na(Pfluoro)]
lines(y=c(min(Pfluoro),max(Pfluoro)), x=c(min(Pfluoro),max(Pfluoro)), lwd=2, col="grey")
summary(fv.asreml.02)
summary(fv.asreml.02[fv.asreml.02>0])
summary(Pfluoro)


data.rewm <- data.frame(stn.id,lat.stn,long.stn,re.stn.wm1,re.stn.wm2,re.stn.wm3,fv.asreml.02.stn,fluoro.stn)

write.csv(x=data.rewm, file="data.rewm.csv")

summary(glm.data$fluoro)

summary((exp(fv.asreml.02)-1))

summary((exp(L.fluoro)-1))


# graph just depth profile

rm(temp1.pv)

# graph profile.depth splines

#depth.v <- min(glm.data$profile.depth)+(max(glm.data$profile.depth)-min(glm.data$profile.depth))*(seq(1,200,1)-1)/199

depth.v <- seq(2,250,2)

temp1.pv <- predict(asreml.02, classify =list("profile.depth"), levels=list("profile.depth"=
     list(profile.depth=depth.v)))

asreml.02.pv<-temp1.pv$predictions$"profile.depth"$pvals

asreml.02.pv.est<-as.double(asreml.02.pv[,2])

asreml.02.pv.seest<-as.double(asreml.02.pv[,3])

summary(asreml.02.pv.est)

summary(asreml.02.pv.seest)

asreml.02.pv.est.m <- t(matrix(data=asreml.02.pv.est, ncol=length(depth.v), nrow=1, byrow=T))

asreml.02.pv.est.mul <- t(matrix(data=(asreml.02.pv.est+2*asreml.02.pv.seest), ncol=length(depth.v), nrow=1, byrow=T))

asreml.02.pv.est.mll <- t(matrix(data=(asreml.02.pv.est-2*asreml.02.pv.seest), ncol=length(depth.v), nrow=1, byrow=T))

asreml.02.pv.est.m <- exp(asreml.02.pv.est.m) - 1

asreml.02.pv.est.mul <- exp(asreml.02.pv.est.mul) - 1

asreml.02.pv.est.mll <- exp(asreml.02.pv.est.mll) - 1

summary(as.vector(asreml.02.pv.est.m))

plot(y=asreml.02.pv.est.m[,1], x=depth.v, ylim=c(-0.2,2.0), xlim=c(0,250), xlab="profile.depth",
    ylab="Predicted Relative Fluoro", type="l", lwd=2, col=1)

lines(y=asreml.02.pv.est.mul[,1], x=depth.v, lty=3, lwd=2, col="black")
lines(y=asreml.02.pv.est.mll[,1], x=depth.v, lty=3, lwd=2, col="black")

plot(y=asreml.02.pv.est.m[,1], x=depth.v, ylim=c(-0.2,25), xlim=c(0,250), xlab="profile.depth",
    ylab="Predicted Relative Fluoro", type="l", lwd=2, col=2)

lines(y=asreml.02.pv.est.mul[,1], x=depth.v, lty=3, lwd=2, col=3)
lines(y=asreml.02.pv.est.mll[,1], x=depth.v, lty=3, lwd=2, col=3)
points(y=glm.data$fluoro, x=glm.data$profile.depth, pch=1, cex=0.4)


# graph just par.ml splines

rm(temp1.pv)

# graph par.ml splines

par.v <- 1*(seq(1,200,1)-1)/199

temp1.pv <- predict(asreml.02, classify =list("profile.depth:par.ml"), levels=list("profile.depth:par.ml"=
     list(profile.depth=40,par.ml=par.v)))

asreml.02.pv<-temp1.pv$predictions$"profile.depth:par.ml"$pvals

asreml.02.pv.est<-as.double(asreml.02.pv[,3])

asreml.02.pv.seest<-as.double(asreml.02.pv[,4])

summary(asreml.02.pv.est)

summary(asreml.02.pv.seest)

asreml.02.pv.est.m <- t(matrix(data=asreml.02.pv.est, ncol=length(par.v), nrow=1, byrow=T))

asreml.02.pv.est.mul <- t(matrix(data=(asreml.02.pv.est+2*asreml.02.pv.seest), ncol=length(par.v), nrow=1, byrow=T))

asreml.02.pv.est.mll <- t(matrix(data=(asreml.02.pv.est-2*asreml.02.pv.seest), ncol=length(par.v), nrow=1, byrow=T))

asreml.02.pv.est.m <- exp(asreml.02.pv.est.m) - 1

asreml.02.pv.est.mul <- exp(asreml.02.pv.est.mul) - 1

asreml.02.pv.est.mll <- exp(asreml.02.pv.est.mll) - 1

summary(as.vector(asreml.02.pv.est.m))

summary(glm.data$par)

plot(y=asreml.02.pv.est.m[,1], x=par.v, ylim=c(1.0,2.0), xlim=c(0,1), xlab="PAR",
    ylab="Predicted Relative Fluoro", type="l", lwd=2, col=1)

lines(y=asreml.02.pv.est.mul[,1], x=par.v, lty=3, lwd=2, col="black")
lines(y=asreml.02.pv.est.mll[,1], x=par.v, lty=3, lwd=2, col="black")

plot(y=asreml.02.pv.est.m[,1], x=par.v, ylim=c(-0.2,2.5), xlim=c(0,1), xlab="PAR",
    ylab="Predicted Relative Fluoro", type="l", lwd=2, col=1)

lines(y=asreml.02.pv.est.mul[,1], x=par.v, lty=3, lwd=2, col="black")
lines(y=asreml.02.pv.est.mll[,1], x=par.v, lty=3, lwd=2, col="black")

points(y=glm.data$fluoro, x=glm.data$par, pch=".")




depth.p.f <- as.factor(glm.data$profile.depth)

par.d.m <- tapply(X=glm.data$par[!is.na(glm.data$par)], INDEX=depth.p.f[!is.na(glm.data$par)], FUN=mean)

par.d.sd <- tapply(X=glm.data$par[!is.na(glm.data$par)], INDEX=depth.p.f[!is.na(glm.data$par)], FUN=var)

par.d.n <- tapply(X=rep(1,sum(!is.na(glm.data$par))), INDEX=depth.p.f[!is.na(glm.data$par)], FUN=sum)

par.UL <- as.vector(par.d.m) + 2*(as.vector(par.d.sd)/as.vector(par.d.n))^0.5
par.LL <- as.vector(par.d.m) - 2*(as.vector(par.d.sd)/as.vector(par.d.n))^0.5


depth.d.m <- tapply(X=glm.data$profile.depth, INDEX=list(depth.p.f), FUN=mean)

summary(par.d.m)
summary(depth.d.m)

plot(y=par.d.m, x=depth.d.m, xlab="depth.profile", ylab="Mean PAR", type="b", ylim=c(0,0.5))
for (i in 1:length(par.d.m)) {
segments(x0=depth.d.m[i], x1=depth.d.m[i], y0=par.LL[i], y1=par.UL[i])}


# graph just temp.ml splines

rm(temp1.pv)

# graph temp.ml splines

temp.v <-  min(glm.data$temp[!is.na(glm.data$temp)])+(max(glm.data$temp[!is.na(glm.data$temp)])
-min(glm.data$temp[!is.na(glm.data$temp)]))*(seq(1,200,1)-1)/199


temp1.pv <- predict(asreml.02, classify =list("profile.depth:temp.ml"), levels=list("profile.depth:temp.ml"=
     list(profile.depth=40,temp.ml=temp.v)))

asreml.02.pv<-temp1.pv$predictions$"profile.depth:temp.ml"$pvals

asreml.02.pv.est<-as.double(asreml.02.pv[,3])

asreml.02.pv.seest<-as.double(asreml.02.pv[,4])

summary(asreml.02.pv.est)

summary(asreml.02.pv.seest)

asreml.02.pv.est.m <- t(matrix(data=asreml.02.pv.est, ncol=length(temp.v), nrow=1, byrow=T))

asreml.02.pv.est.mul <- t(matrix(data=(asreml.02.pv.est+2*asreml.02.pv.seest), ncol=length(temp.v), nrow=1, byrow=T))

asreml.02.pv.est.mll <- t(matrix(data=(asreml.02.pv.est-2*asreml.02.pv.seest), ncol=length(temp.v), nrow=1, byrow=T))

asreml.02.pv.est.m <- exp(asreml.02.pv.est.m) - 1

asreml.02.pv.est.mul <- exp(asreml.02.pv.est.mul) - 1

asreml.02.pv.est.mll <- exp(asreml.02.pv.est.mll) - 1

summary(as.vector(asreml.02.pv.est.m))

summary(glm.data$temp)

plot(y=asreml.02.pv.est.m[,1], x=temp.v, ylim=c(0.0,2.0), xlim=c(-2,10), xlab="Temperature",
    ylab="Predicted Relative Fluoro", type="l", lwd=2, col=1)

lines(y=asreml.02.pv.est.mul[,1], x=temp.v, lty=3, lwd=2, col="black")
lines(y=asreml.02.pv.est.mll[,1], x=temp.v, lty=3, lwd=2, col="black")



# graph just ice.fdays.ml splines

rm(temp1.pv)

# graph ice.free.days.ml splines

ice.fdays <- glm.data$ice.fdays

ice.fdays.v <- min(ice.fdays)+(max(ice.fdays)-min(ice.fdays))*(seq(1,200,1)-1)/199


temp1.pv <- predict(asreml.02, classify =list("profile.depth:ice.fdays.ml"), levels=list("profile.depth:ice.fdays.ml"=
     list(profile.depth=40,ice.fdays.ml=ice.fdays.v)))

asreml.02.pv<-temp1.pv$predictions$"profile.depth:ice.fdays.ml"$pvals

asreml.02.pv.est<-as.double(asreml.02.pv[,3])

asreml.02.pv.seest<-as.double(asreml.02.pv[,4])

summary(asreml.02.pv.est)

summary(asreml.02.pv.seest)

asreml.02.pv.est.m <- t(matrix(data=asreml.02.pv.est, ncol=length(ice.fdays.v), nrow=1, byrow=T))

asreml.02.pv.est.mul <- t(matrix(data=(asreml.02.pv.est+2*asreml.02.pv.seest), ncol=length(ice.fdays.v), nrow=1, byrow=T))

asreml.02.pv.est.mll <- t(matrix(data=(asreml.02.pv.est-2*asreml.02.pv.seest), ncol=length(ice.fdays.v), nrow=1, byrow=T))

asreml.02.pv.est.m <- exp(asreml.02.pv.est.m) - 1

asreml.02.pv.est.mul <- exp(asreml.02.pv.est.mul) - 1

asreml.02.pv.est.mll <- exp(asreml.02.pv.est.mll) - 1

summary(as.vector(asreml.02.pv.est.m))

summary(glm.data$par)

plot(y=asreml.02.pv.est.m[,1], x=ice.fdays.v, ylim=c(0.0,5.0), xlab="Ice free days",
    ylab="Predicted Relative Fluoro", type="l", lwd=2, col=1)

lines(y=asreml.02.pv.est.mul[,1], x=ice.fdays.v, lty=3, lwd=2, col="black")
lines(y=asreml.02.pv.est.mll[,1], x=ice.fdays.v, lty=3, lwd=2, col="black")

plot(y=asreml.02.pv.est.m[,1], x=ice.fdays.v, ylim=c(-0.2,2.5), xlim=c(0,1), xlab="PAR",
    ylab="Predicted Relative Fluoro", type="l", lwd=2, col=1)

lines(y=asreml.02.pv.est.mul[,1], x=ice.fdays.v, lty=3, lwd=2, col="black")
lines(y=asreml.02.pv.est.mll[,1], x=ice.fdays.v, lty=3, lwd=2, col="black")

points(y=glm.data$fluoro, x=glm.data$ice.fdays, pch=".")



# graph water.mass:ice.fdays splines

ice.fdays <- glm.data$ice.fdays

summary(ice.fdays)

ice.fdays.v <- min(ice.fdays)+(max(ice.fdays)-min(ice.fdays))*(seq(1,200,1)-1)/199

rm(temp1.pv)

temp1.pv <- predict(asreml.02, classify = list("profile.depth:water.mass.f:ice.fdays.ml"),
        levels=list("profile.depth:water.mass.f:ice.fdays.ml"= list(profile.depth=40,water.mass.f=levels(glm.data$water.mass.f),ice.fdays.ml=ice.fdays.v)))

length(temp1.pv$predictions)

asreml.02.pv<-temp1.pv$predictions$"profile.depth:water.mass.f:ice.fdays.ml"$pvals

dim(asreml.02.pv)

asreml.02.pv.est<-as.double(asreml.02.pv[,4])

asreml.02.pv.seest<-as.double(asreml.02.pv[,5])

summary(asreml.02.pv.est)

summary(asreml.02.pv.seest)

length(ice.fdays.v)

asreml.02.pv.est.m <- t(matrix(data=asreml.02.pv.est, ncol=length(ice.fdays.v), nrow=4, byrow=F))

asreml.02.pv.seest.m <- t(matrix(data=asreml.02.pv.seest, ncol=length(ice.fdays.v), nrow=4, byrow=F))

asreml.02.pv.est.mul <- t(matrix(data=(asreml.02.pv.est+2*asreml.02.pv.seest), ncol=length(ice.fdays.v), nrow=4, byrow=F))

asreml.02.pv.est.mll <- t(matrix(data=(asreml.02.pv.est-2*asreml.02.pv.seest), ncol=length(ice.fdays.v), nrow=4, byrow=F))

asreml.02.pv[1:200,]

asreml.02.pv.est.m[1:200,1]


asreml.02.pv.est.m <- exp(asreml.02.pv.est.m) - 1

asreml.02.pv.est.mul <- exp(asreml.02.pv.est.mul) - 1

asreml.02.pv.est.mll <- exp(asreml.02.pv.est.mll) - 1

dim(asreml.02.pv.est.m)

summary(as.vector(asreml.02.pv.est.m))

summary(c(as.vector(asreml.02.pv.est.mul),as.vector(asreml.02.pv.est.mll)))

plot(y=asreml.02.pv.est.m[,2], x=ice.fdays.v, ylim=c(0.0,5.0), xlim=c(-35,105), xlab="ice free days",
    ylab="Predicted Relative Fluoro", type="l", lwd=2)

lines(y=asreml.02.pv.est.mul[,2], x=ice.fdays.v, lty=2, col=1, lwd=2)
lines(y=asreml.02.pv.est.mll[,2], x=ice.fdays.v, lty=2, col=1, lwd=2)

lines(y=asreml.02.pv.est.m[,2], x=ice.fdays.v, lty=3, col=2)
lines(y=asreml.02.pv.est.m[,3], x=ice.fdays.v, lty=4, col=3)
lines(y=asreml.02.pv.est.m[,4], x=ice.fdays.v, lty=5, col=4)

legend(y=2.0, x=0.6, legend=c("Wm1","Wm2","Wm3","Wm4","Wm1.95%cl"), lty=c(1,3,4,5,2), col=c(1,2,3,4,1))


# graph par splines

par.v <- min(par.m)+(max(par.m)-min(par.m))*(seq(1,200,1)-1)/199

rm(temp1.pv)

temp1.pv <- predict(asreml.02, classify = list("water.mass:ice.free.days.ml"),
        levels=list("water.mass:ice.free.days.ml"= list(water.mass=levels(water.mass),ice.free.days.ml=par.v)))

length(temp1.pv$predictions)

asreml.02.pv<-temp1.pv$predictions$"water.mass:ice.free.days.ml"$pvals

dim(asreml.02.pv)

asreml.02.pv.est<-as.double(asreml.02.pv[,3])

asreml.02.pv.seest<-as.double(asreml.02.pv[,4])

summary(asreml.02.pv.est)

summary(asreml.02.pv.seest)

length(ice.free.days.ml)

asreml.02.pv.est.m <- t(matrix(data=asreml.02.pv.est, ncol=length(par.v), nrow=4, byrow=F))

asreml.02.pv.seest.m <- t(matrix(data=asreml.02.pv.seest, ncol=length(par.v), nrow=4, byrow=F))

asreml.02.pv.est.mul <- t(matrix(data=(asreml.02.pv.est+2*asreml.02.pv.seest), ncol=length(par.v), nrow=4, byrow=F))

asreml.02.pv.est.mll <- t(matrix(data=(asreml.02.pv.est-2*asreml.02.pv.seest), ncol=length(par.v), nrow=4, byrow=F))

asreml.02.pv[1:20,]

asreml.02.pv.est.m[1:20,1]


asreml.02.pv.est.m <- exp(asreml.02.pv.est.m) - 1

asreml.02.pv.est.mul <- exp(asreml.02.pv.est.mul) - 1

asreml.02.pv.est.mll <- exp(asreml.02.pv.est.mll) - 1



summary(as.vector(asreml.02.pv.est.m))

summary(c(as.vector(asreml.02.pv.est.mul),as.vector(asreml.02.pv.est.mll)))

graphsheet( )

plot(y=asreml.02.pv.est.m[,1], x=par.v, ylim=c(-0.6,1.5), xlim=c(0,1.0), xlab="par",
    ylab="Predicted Relative Fluoro", type="l")

lines(y=asreml.02.pv.est.mul[,1], x=par.v, lty=2, col=1)
lines(y=asreml.02.pv.est.mll[,1], x=par.v, lty=2, col=1)

lines(y=asreml.02.pv.est.m[,2], x=par.v, lty=3, col=2)
lines(y=asreml.02.pv.est.m[,3], x=par.v, lty=4, col=3)

legend(y=1.5, x=0.6, legend=c("Wm1","Wm2","Wm3","Wm1.95%cl"), lty=c(1,3,4,2), col=c(1,2,3,1))

lines(y=asreml.02.pv.est.m[,4], x=par.v, lty=5, col=4)

legend(y=1.5, x=0.6, legend=c("Wm1","Wm2","Wm3","Wm4","Wm1.95%cl"), lty=c(1,3,4,5,2), col=c(1,2,3,4,1))



