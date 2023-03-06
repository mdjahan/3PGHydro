library(devtools)
install_github("mdjahan/3PGHydro/rpackage.3PGHydro/")
library(rpackage.3PGHydro)
?run_3PGhydro
#
setwd("C:/Users/Marc Djahangard/Documents/IFE/2022/3_3PG/3_3PG_Hydro_R/example")

#climate data
climate <- read.csv("climate.csv")
climate$date <- as.Date(climate$date,format="%d/%m/%Y")

#load 3PG parameters
Parameter <- read.csv("3PG_Parameter.csv") #parameters from Forrester et al. 2021
#select species
p <- Parameter$Quercus.petraea

#Site & Stand characteristics
lat <- 47.5
StartDate <- "01/01/1960"
StandAgei <- 20
EndAge <- 65
StemNoi <- 500
WSi <- 200
WFi <- 10
WRi <- 20
CO2Concentration <- "Historical"
FR <- 0.5
HeightEquation <- 1 
SVEquation <- 2
#Soil
SoilClass <- 2
EffectiveRootZoneDepth <- 1
DeepRootZoneDepth <- 5
RocksER <- 0.2 
RocksDR <- 0.4 
thinAges <- NULL #c(30,40,50,60)
thinVals <- NULL #c(450,400,350,300)
thinWF <- NULL #c(1,1,1,1)
thinWR <- NULL #c(1,1,1,1)
thinWS <- NULL #c(1,1,1,1)

OutputRes <- "daily"

#Run 3PG-Hydro
out <- run_3PGhydro(climate,p,lat,StartDate,StandAgei,EndAge,WFi,WRi,WSi,StemNoi,CO2Concentration,FR,SVEquation,HeightEquation,SoilClass,EffectiveRootZoneDepth,DeepRootZoneDepth,RocksER,RocksDR,thinAges,thinVals,thinWF,thinWR,thinWS,OutputRes)

#Yearly Output
OutputRes <- "yearly"
out_yearly <-  run_3PGhydro(climate,p,lat,StartDate,StandAgei,EndAge,WFi,WRi,WSi,StemNoi,CO2Concentration,FR,SVEquation,HeightEquation,SoilClass,EffectiveRootZoneDepth,DeepRootZoneDepth,RocksER,RocksDR,thinAges,thinVals,thinWF,thinWR,thinWS,OutputRes)

