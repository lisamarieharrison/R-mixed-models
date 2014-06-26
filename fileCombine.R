wd='C:\\Users\\Lisa\\Documents\\phd\\southern ocean\\BROKE-West_CTD_data\\Processed\\CTDData_200506031_080704'
maxRow=125
outFn='procCTD.csv'
oxygenDontCareVal=-9

files=list.files(wd,pattern='.csv')
dat=read.csv(paste(wd,files[1],sep='\\'))
dat=dat[1:maxRow,]
write.csv(dat,paste(wd,outFn,sep='\\'),row.names=FALSE)

for(i in 2:length(files))
{
  tmp=read.csv(paste(wd,files[i],sep='\\'))
  tmp=tmp[1:maxRow,]
  message(i,' ; ', ncol(tmp),'\n')
  write.table(tmp,paste(wd,outFn,sep='\\'),row.names=FALSE,col.names=FALSE,append=TRUE,sep=',') 
}

dat=read.csv(paste(wd,outFn,sep='\\'))
dat=dat[-which(is.na(dat$Cast.Number)),]
dat$Oxygen..micromole.litre..micromole.litre[dat$Oxygen..micromole.litre..micromole.litre==oxygenDontCareVal]=NA
names(dat)
table(dat$Cast.Number)

subdat=subset(dat,Cast.Number==1)
plot(subdat$Oxygen..micromole.litre..micromole.litre,subdat$Pressure,ylim=c(250,0),
     xlim=range(dat$Oxygen..micromole.litre..micromole.litre,na.rm=TRUE),type='o')
for(i in 2:max(dat$Cast.Number)){
  subdat=subset(dat,Cast.Number==i)
  points(subdat$Oxygen..micromole.litre..micromole.litre,subdat$Pressure,col='grey')
}

write.table(dat,paste(wd,outFn,sep='\\'),row.names=FALSE,col.names=TRUE,append=FALSE,sep=',') 
