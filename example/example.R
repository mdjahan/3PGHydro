#library(devtools)
#install_github("mdjahan/3PGHydro/rpackage.3PGHydro/")
library(rpackage.3PGHydro)
?run_3PGhydro
#
setwd("YOURPATH")

#climate data
climate <- read.csv("climate.csv")
climate$date <- as.Date(climate$date,format="%d/%m/%Y")

#load 3PG parameters
Parameter <- read.csv("parameter.csv") #P.abies from Forrester et al. 2021; F. sylvatica from Augustynczik et al 2017
#select species
p <- Parameter$Fagus.sylvatica

#Site & Stand characteristics
lat <- 47.5
StartDate <- "01/01/1960"
StandAgei <- 30
EndAge <- 80
StemNoi <- 600
WSi <- 135
WFi <- 8
WRi <- 15
CO2Concentration <- "Historical"
FR <- 0.5
HeightEquation <- 1 
SVEquation <- 1
#Soil
SoilClass <- 2
EffectiveRootZoneDepth <- 1
DeepRootZoneDepth <- 3
RocksER <- 0.2 
RocksDR <- 0.4 

#Management
thinAges <- NULL #c(30,40,50,60)
thinVals <- NULL #c(450,400,350,300)
thinWF <- NULL #c(1,1,1,1)
thinWR <- NULL #c(1,1,1,1)
thinWS <- NULL #c(1,1,1,1)

#Output Resolution
OutputRes <- "daily"

#RUN: 3PG-Hydro
out <- run_3PGhydro(climate,p,lat,StartDate,StandAgei,EndAge,WFi,WRi,WSi,StemNoi,CO2Concentration,FR,HeightEquation,SVEquation,SoilClass,EffectiveRootZoneDepth,DeepRootZoneDepth,RocksER,RocksDR,thinAges,thinVals,thinWF,thinWR,thinWS,OutputRes)

#Some Plots:
par(mfrow=c(2,2))
plot(y=out$Height,x=out$Date,main="Height")
plot(y=out$avDBH,x=out$Date,main="DBH")
plot(y=out$StandVol,x=out$Date,main="Standing Volume")
plot(y=out$StemNo,x=out$Date,main="Stems")

par(mfrow=c(2,2))
plot(out$GPP[730:1095],x=out$Date[730:1095],type="l",main="GPP")
plot(out$volWCer[730:1095],x=out$Date[730:1095],type="l",main="vol WC ER")
plot(out$DeepPercolation[730:1095],x=out$Date[730:1095],type="l",main="Deep Percolation")
plot(out$LAI[730:1095],x=out$Date[730:1095],type="l",main="LAI")


#Yearly Output
OutputRes <- "yearly"
out_yearly <-  run_3PGhydro(climate,p,lat,StartDate,StandAgei,EndAge,WFi,WRi,WSi,StemNoi,CO2Concentration,FR,HeightEquation,SVEquation,SoilClass,EffectiveRootZoneDepth,DeepRootZoneDepth,RocksER,RocksDR,thinAges,thinVals,thinWF,thinWR,thinWS,OutputRes)

