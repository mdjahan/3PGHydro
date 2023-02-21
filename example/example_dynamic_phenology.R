library(devtools)
install_github("mdjahan/3PGHydro/rpackage.3PGHydro/")
library(rpackage.3PGHydro)
?run_3PGhydro
#
setwd("C:/Users/anoelte/Desktop/Sammelsurium/3PG/3PGHydro/example/")

#climate data
climate <- read.csv("climate_kl.csv")
#climate$date <- as.Date(climate$date,"%d/%m/%Y")

#load 3PG parameters
Parameter <- read.csv("Parameter_Quercus.csv") #parameters from Forrester et al. 2021 and Noelte et al. 2020 and Marc

#select parameter set
Parameter$Noelte.dynamicPhenology[16] <- 1.38
Parameter$Noelte.dynamicPhenology[17] <- 0.7

Parameter$Forrester[1] <- 1.12
Parameter$Forrester[2] <- 0.03
Parameter$Forrester[15] <- 18
Parameter$Forrester[16] <- 1.29
Parameter$Forrester[17] <- 0.72
Parameter$Forrester[40] <- 0.042

p <- Parameter$Noelte.dynamicPhenology
p <- Parameter$Forrester
p <- Parameter$Marc

#Site & Stand characteristics
lat <- 49.39
StartDate <- "01/11/1990"
StandAgei <- 54
EndAge <- 79
StemNoi <- 1880
WSi <- 122
WFi <- 5
WRi <- 22
CO2Concentration <- "Historical"
FR <- 0.6
HeightEquation <- 1
SVEquation <- 2 #1 is with parameters and form SV = a * H^b * D^c (Forrester et al. 2021), 2: with WS & density, 3: backup --> simple with form factor 0.7
#Soil
SoilClass <- 1
EffectiveRootZoneDepth <- 1
DeepRootZoneDepth <- 5
RocksER <- 0.2 #share of rocky parts/skeleton in effective root zone
RocksDR <- 0.4 #share of rocky parts/skeleton in deep root zone
thinAges <- c(59.02, 64.02, 69.02, 74.02, 79.02) #thinning after inventory
thinVals <- c(1156,  780,  656,  556,  484) 
thinWF <- c(0.7, 0.9, 0.4, 0.4, 0.2) 
thinWR <- c(0.7, 0.9, 0.4, 0.4, 0.2)  
thinWS <- c(0.7, 0.9, 0.4, 0.4, 0.2)
phenology <- 0

#Run the function
#source('C:/Users/anoelte/Desktop/Sammelsurium/3PG/3PGHydro/rpackage.3PGHydro/R/run_3PGhydro.R')
out <- run_3PGhydro(climate,p,lat,StartDate,StandAgei,EndAge,WFi,WRi,WSi,StemNoi,CO2Concentration,FR,SVEquation,SoilClass,EffectiveRootZoneDepth,DeepRootZoneDepth,RocksER,RocksDR,thinAges,thinVals,thinWF,thinWR,thinWS,phenology)


#Plot
library(ggplot2)
ggplot()+
  geom_line(aes(x=StandAge, y = WS), data = out)

library(dplyr)
out$year <- as.numeric(format(out$Date,"%Y"))
npp_yearly <- out %>% group_by(year) %>% summarize(npp = sum(NPP, na.rm = T))
#with both paramters sets (forrester and nölte), npp is much too low:
  # with approx. 3-5 t/ha/yr, but it should be between 7-12 t/ha/yr
  # about 4 at the beginning and 1 at the end of the simulation period 
  # since NPP is low not enough Carbon for WFstore and WF decreses in a negative feedback loop
  # with the Forrester parameters almost no leaf biomass is allocated to WF, almost everything to stem (very small pFS20 parameter value)
  # this is why WF remains stable and stem biomass increase is still reasonable, albeit too low, with the low NPP production
  # WF is kind of independent of everything, except for the early years (high pFS2 parameter value)

#Test with no thinning
thinAges <- NULL
thinVals <- NULL
thinWF <- NULL 
thinWR <- NULL 
thinWS <- NULL 
