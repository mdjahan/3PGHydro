library(devtools)
install_github("mdjahan/3PGHydro/rpackage.3PGHydro/")
library(rpackage.3PGHydro)
?run_3PGhydro
#
setwd("YOUR_DATA_PATH")

#climate data
climate <- read.csv("climate.csv")

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
HeightEquation <- 2 
SVEquation <- 1 #1 is with parameters and form SV = a * H^b * D^c (Forrester et al. 2021), 2: with WS & density, 3: backup --> simple with form factor 0.7
#Soil
SoilClass <- 2
EffectiveRootZoneDepth <- 1
DeepRootZoneDepth <- 5
RocksER <- 0.2 #share of rocky parts/skeleton in effective root zone
RocksDR <- 0.4 #share of rocky parts/skeleton in deep root zone
thinAges <- NULL
thinVals <- NULL
thinWF <- NULL 
thinWR <- NULL 
thinWS <- NULL 

#Run the function
out <- run_3PGhydro(climate,p,lat,StartDate,StandAgei,EndAge,WFi,WRi,WSi,StemNoi,CO2Concentration,FR,SVEquation,SoilClass,EffectiveRootZoneDepth,DeepRootZoneDepth,RocksER,RocksDR,thinAges,thinVals,thinWF,thinWR,thinWS)



