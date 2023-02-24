#' Run 3PG-Hydro
#'
#' 3PG-Hydro is an update of the original 3PG Forest Growth model by Landsberg and Waring (1997) (DOI:10.1016/S0378-1127(97)00026-1). 3PG-Hydro calculates on a daily timestep, includes a soil-water-model as well as a snow routine. Further information in the publication on 3PG-Hydro by Yousefpour and Djahangard (2021) (DOI:10.3390/f12121729). 3PG-Hydro is available as an R-package, coding done by Anja Nölte & Marc Djahangard.
#' @param climate Climate data as .csv file, mandatory column names ("date" ("dd/mm/yyyy"),"Tav","Tmax","Tmin","Rain","SolarRad")
#' @param p 3PG tree species parameter .csv file
#' @param lat site latitude UTM
#' @param StartDate starting date, format: "dd/mm/yyyy"
#' @param StandAgei starting stand age in years
#' @param EndAge end stand age in years
#' @param WFi inital weight of foliage in tons
#' @param WRi inital weight of roots in tons
#' @param WSi inital weight of stems in tons
#' @param StemNoi initial stem number
#' @param CO2Concentration select which CO2 concentration equation should be used: "Historical", "RCP2.6", "RCP4.5", "RCP8.5"
#' @param FR fertility rating of site (0-1)
#' @param HeightEquation choose the height equation (1-2), (1: 3PG Original Equation, 2: Michajlow-Schumacher from Forrester et al. (2021))
#' @param SVEquation choose the equation used for stand volume calculation (1-3) (1: V = aV * H^VH * D^VB, Forrester et al. (2021), 2: 3PG original function with stem mass & density, 3: simple volume calculation V = FormFactor * CylinderVolume)
#' @param SoilClass soil classes, type number 1 - 4 (1: sand, 2: sandy loam, 3: clay loam, 4: clay)
#' @param EffectiveRootZoneDepth depth of effective root zone in meter
#' @param DeepRootZoneDepth depth of deep root zone in meter
#' @param RocksER share of rocks/skeleton in effective root zone (0-1) 0:none, 1:all
#' @param RocksDR share of rocks/skeleton in deep root zone (0-1) 0:none, 1:all
#' @param thinAges ages at which thinning should be applied, in case of no thinning: NULL
#' @param thinVals number of stems left after thinning, in case of no thinning: NULL
#' @param thinWF thinning regime foliage below to above (0.5 - 1.5), in case of no thinning: NULL
#' @param thinWR thinning regime foliage below to above (0.5 - 1.5), in case of no thinning: NULL
#' @param thinWS thinning regime foliage below to above (0.5 - 1.5), in case of no thinning: NULL
#' @return output file in daily time steps: Date, StandAge, StemNo, WF, WR, WS, DBH, Height, StandVol, volWCer, volWCdr, NPP, NEE, LAI, Evapotranspiration, AvStemMass, Basal Area, Self Thinning, WSext, StandVol_loss, VolProduction_tot, Deep Percolation, Run Off 
#' @examples 
#'climate <- read.csv("climate.csv")
#'p <- read.csv("3PG_Parameter.csv") 
#'p <- p$species1 #select species parameters column
#' lat <- 40.5
#' StartDate <- "01/01/1980"
#' StandAgei <- 40
#' EndAge <- 65
#' StemNoi <- 300
#' WSi <- 150
#' WFi <- 10
#' WRi <- 20
#' CO2Concentration <- "Historical"
#' FR <- 0.7
#' HeightEquation <- 1
#' SVEquation <- 2
#' SoilClass <- 2
#' EffectiveRootZoneDepth <- 1
#' DeepRootZoneDepth <- 5
#' RocksER <- 0.2 
#' RocksDR <- 0.4 
#' thinAges <- c(50,55,60)
#' thinVals <- c(250,200,150) 
#' thinWF <- rep(0.8,length(thinAges)) 
#' thinWR <- rep(0.8,length(thinAges)) 
#' thinWS <- rep(0.8,length(thinAges))  
#' out <- run_3PGhydro(climate,p,lat,StartDate,StandAgei,EndAge,WFi,WRi,WSi,StemNoi,CO2Concentration,FR,SoilClass,EffectiveRootZoneDepth,DeepRootZoneDepth,RocksER,RocksDR,thinAges,thinVals,thinWF,thinWR,thinWS)
#' @export
run_3PGhydro <- function(climate,p,lat,StartDate,StandAgei,EndAge,WFi,WRi,WSi,StemNoi,CO2Concentration,FR,HeightEquation,SVEquation,SoilClass,EffectiveRootZoneDepth,DeepRootZoneDepth,RocksER,RocksDR,thinAges,thinVals,thinWF,thinWR,thinWS){
  
  ############################################################
  #parameters
  ##############################################################
  pFS2 <- p[1]
  pFS20 <- p[2]
  aWs <- p[3]
  nWs <- p[4]
  pRx <- p[5]
  pRn <- p[6]
  gammaF1 <- p[7]
  gammaF0 <- p[8]
  tgammaF <- p[9]
  gammaR <- p[10]
  leafgrow <- p[11]
  leaffall <- p[12]
  Tmin <- p[13]
  Tmax <- p[14]
  Topt <- p[15]
  fCalpha700 <- p[16]
  fCg700 <- p[17]
  m0 <- p[18]
  fN0 <- p[19]
  fNn <- p[20]
  MaxAge <- p[21]
  nAge <- p[22]
  rAge <- p[23]
  gammaN1 <- p[24]
  gammaN0 <- p[25]
  tgammaN <- p[26]
  ngammaN <- p[27]
  wSx1000 <- p[28]
  thinPower <- p[29]
  mF <- p[30]
  mR <- p[31]
  mS <- p[32]
  SLA0 <- p[33]
  SLA1 <- p[34]
  tSLA <- p[35]
  k <- p[36]
  fullCanAge <- p[37]
  MaxIntcptn <- p[38]
  LAImaxIntcptn <- p[39]
  alphaCx <- p[40]
  y <- p[41]
  CoeffCond <- p[42] #defines stomatal response to VPD
  MinCond <- p[43]
  MaxCond <- p[44] #maximum stomatal conductance integrated over the canopy
  BLcond <- p[45] #conductance between the canopy layer and the atmoshperic layer
  LAIgcx <- p[46] #canopy LAI for maximum canopy conductance (stomata conductance)
  fracBB0 <- p[47]
  fracBB1 <- p[48]
  tBB <- p[49]
  rho0 <- p[50]
  rho1 <- p[51]
  tRho <- p[52]
  aH <- p[53]
  nHB <- p[54]
  nHC <- p[55]
  aV <- p[56]
  nVB <- p[57]
  nVH <- p[58]
  gammaN0attack <- p[59]
  attackAge <- p[60]
  attackTime <- p[61] 
  #
  kF <- 1
  SnowmMeltFactor <- 2.5
  #Conversion factors
  gDM_mol <- 24 #conversion of mol to gDM (Molecular weight of dry matter)
  molPAR_MJ <- 2.3 #conversion of MJ to PAR (mol/m2)
  Qa <- -90 #Intercept of net v. solar radiation relationship
  Qb <- 0.8 # slope o net v. solar radiation relationship
  #Parameters Penman-Monteith equation for computing canopy transpiration
  #The following are constants in the PM formula (Landsberg & Gower, 1997)
  e20 <- 2.2          #rate of change of saturated VP with T at 20C
  rhoAir <- 1.2       # density of air, kg/m3
  lambda <- 2460000#  # latent heat of vapourisation of H2O (J/kg)
  VPDconv <- 0.000622 # convert VPD to saturation deficit = 18/29/1000
  #Parameter for Soil Evaporation
  maxgSoil <- 0.001 #maximum soil conductance
  gAS <- 0.015 #Soil aerodynamic conductance
  #For NEE calculation (from Meyer et al. (2018)
  Rh10soil <- 0.00025 #base-respiration rate at 10 degrees (tonsDM/tonC/d); range: 0.0001 - 0.0005
  Q10 <- 2 #Q10 for heterotrophic respiration
  Csoil <- 40 #soil carbon content (tons/ha)
  #  
  poolFractn <- 0
  
  #CO2 Equations: fitted with values from IIASA
  CO2HistEq <- function(year) { #From 1850 - 2020
    CO2 <- -3.7892295404645533e+005+5.9852877232176991e+002*year+(-3.1497736226299061e-001)*year^2+5.5268027924225127e-005*year^3
    return(CO2)
  }
  CO2rcp26Eq <- function(year) {
    CO2 <- -0.014632*year^2+60.2673*year+61617
    return(CO2)
  }
  CO2rcp45Eq <- function(year) {
    CO2 <- -0.0200214*year^2+84.1081*year-87793.8
    return(CO2)
  }
  CO2rcp85Eq <- function(year) {
    CO2 <- 0.044329*year^2-176.052*year+175158
    return(CO2)
  }
  
  #####################################################
  #Initialisation
  ####################################################
  
  #Assign initial state of stand
  StandAge <-StandAgei
  EndAge <- EndAge
  StartDate <- as.Date(StartDate,"%d/%m/%Y")
  EndDate <-   seq(StartDate, length=(EndAge-StandAgei)+1, by="years")[(EndAge-StandAgei)+1]
  date <- StartDate
  year <- as.numeric(format(date,"%Y"))
  currentMonth <- as.numeric(format(as.Date(date,format="%d/%m/%Y"),"%m"))
  Duration <- as.numeric(EndDate-StartDate) 
  WS <- WSi
  WF <- WFi
  WR <- WRi
  TotalW <- WSi + WFi + WRi
  StemNo <- StemNoi
  thinEventNo <- 1
  WSext <- 0
  StandVol_loss <- 0
  WSselfThin <- 0
  WSMort <- 0
  accV <- 0
  #StandVol_ext_single <- 0
  
  #Assign Soil Parameters and 3PG original SWconstant and SWpower as function of Soil Class
  #Soil Class: 1 = sand, 2 = sandy loam, 3 = clay loam, 4 = clay
  if (SoilClass > 0) {
    #3PG Hydro Soil Class Parameters
    if (SoilClass == 1){
      SWconst	<- 0.8
      SWpower	<- 12
      volRes <- 0.02
      volSat <- 0.38
      VGn <- 1.55
      VGalpha <- 4
      kS <- 3.5
      fKdr <- 1
      maxInf <- 30
    }
    if (SoilClass == 2){
      SWconst	<- 0.7
      SWpower	<- 9
      volRes <- 0.08
      volSat <- 0.4
      VGn <- 1.35
      VGalpha <- 3.5
      kS <- 1
      fKdr <- 0.9
      maxInf <- 25
    }
    if (SoilClass == 3){
      SWconst	<- 0.5
      SWpower	<- 5
      volRes <- 0.1
      volSat <- 0.44
      VGn <- 1.25
      VGalpha <- 2.8
      kS <- 0.4
      fKdr <- 0.7
      maxInf <- 20
    }
    if (SoilClass == 4){
      SWconst	<- 4
      SWpower	<- 0.4
      volRes <- 0.12
      volSat <- 0.5
      VGn <- 1.1
      VGalpha <- 2.4
      kS <- 0.1
      fKdr <- 0.5
      maxInf <- 15
    }
  }
  #Soil Depth
  depthER <- EffectiveRootZoneDepth
  depthDR <- DeepRootZoneDepth
  SkER <- RocksER
  SkDR <- RocksDR
  #Soil Wilting Point & Field Capacity
  #h in [m] at FC: 3.36514; at WP: 152.961
  volFC <- (1 / (1 + (VGalpha * 3.36514) ^ VGn)) ^ (1 - 1 / VGn)
  volFC <- volFC * (volSat - volRes) + volRes
  volWP <- (1 / (1 + (VGalpha * 152.961) ^ VGn)) ^ (1 - 1 / VGn)
  volWP <- volWP * (volSat - volRes) + volRes
  #Convert Water Content from volumteric into mm
  #Effective Root Zone:
  tvER <- depthER * 1000 #total volume (liter/m2)
  tvER <- tvER * (1 - SkER) #Subtract Skeleton Share
  erASWsat <- volSat * tvER
  erASWfc <- volFC * tvER
  erASWwp <- volWP * tvER 
  erASWres <- volRes * tvER
  #Deep Root Zone:
  tvDR <- (depthDR - depthER) * 1000
  tvDR <- tvDR * (1 - SkDR)
  drASWsat <- volSat * tvDR #potential add a factor to reduce saturation due to compaction, e.g. 0.9
  drASWfc <- volFC * tvDR #potential add a factor to reduce field capacity due to compaction, e.g. 0.9
  kSdr <- kS * fKdr  
  drASWwp <- volWP * tvDR
  
  #Initial ASW: inital conditions is field capacity
  erASW <- erASWfc
  drASW <- drASWfc
  if(erASW > erASWsat) erASW <- erASWsat
  if(drASW > drASWsat) drASW <- drASWsat
  volWer <- erASW/tvER
  volWdr <- drASW/tvDR
  
  #check fN(FR) for no effect: fNn = 0 ==> fN(FR)=1 for all FR
  if (fNn == 0) fN0 <- 1
  
  #Allocation parameters
  pfsPower <-  log(pFS20 / pFS2) / log(20 / 2)
  pfsConst <- pFS2 / 2 ^ pfsPower
  
  #CO2 Parameters
  fCalphax <- fCalpha700 / (2 - fCalpha700)
  fCg0 <- fCg700 / (2 * fCg700 - 1)
  
  #Initialize soil conditions
  poolFractn <- 0
  applIrrig <- 0
  pooledSW <- 0
  
  #Snow conditions
  aSnow <- 0
  aSnowInt <- 0
  
  #set age depndent factores
  
  #SLA = expF(StandAge, SLA0, SLA1, tSLA, 2)
  if (tSLA != 0) {
    SLA <- SLA1 + (SLA0 - SLA1)*exp(-log(2)*(StandAge/tSLA)^2) 
  } else {
    SLA <- SLA1
  }
  
  #branch and bark fraction
  if (tBB != 0) {
    fracBB <- fracBB1 + (fracBB0 - fracBB1)*exp(-log(2)*(StandAge/tBB)) 
  } else {
    fracBB <- fracBB1
  }
  
  #Density = expF(StandAge, rho0, rho1, tRho, 1)
  if (tRho != 0) {
    Density <- rho1 + (rho0 - rho1)*exp(-log(2)*(StandAge/tRho)) 
  } else {
    Density <- rho1
  }
  
  #gammaF <- gammaFoliage(StandAge)
  if (tgammaF * gammaF1 == 0) {
    gammaF <- gammaF1 
  } else {
    kgammaF <- 12 * log(1 + gammaF1 / gammaF0) / tgammaF #recheck!
    gammaF <- gammaF1 * gammaF0 / (gammaF0 + (gammaF1 - gammaF0) * exp(-kgammaF * StandAge))
  }
  #Leaf fall
  if(leaffall>0){
    currentMonth <- as.numeric(format(as.Date(date,format="%d/%m/%Y"),"%m"))
    if(leafgrow > currentMonth | currentMonth >= leaffall){
      WFprior <- WF
      WF <- 0
    }
  }
  
  #Initialize stand data
  AvStemMass <- WS * 1000 / StemNo  
  avDBH <- (AvStemMass / aWs) ^ (1 / nWs) 
  BasArea <- (((avDBH / 200) ^ 2) * round(pi,9)) * StemNo #Correct for rounding difference of pi
  if(HeightEquation==1) Height <-  aH * avDBH ^ nHB * StemNo ^ nHC
  if(HeightEquation==2) Height <- 1.3 + aH * exp(-nHB/avDBH) + nHC * Density * avDBH #Michajlow-Schumacher (Forrester et. al, 2021)
  LAI <- 0.1*SLA*WF
  #Stand volume m3/ha (excluding branch and bark fraction, see. sands 2002 p.5)
  if(SVEquation==1) StandVol <- aV * (avDBH ^ nVB) * (Height ^ nVH) * StemNo #equation with parameters after Forrester et al. 2021
  if(SVEquation==2) StandVol <- WS * (1 - fracBB) / Density #3PG original equation
  if(SVEquation==3) StandVol <- 0.7 * ((avDBH/200)^2*pi) * Height * StemNo #simple volume equation
  oldV <- StandVol
  
  #Write first line of output (= Start conditions)
  out <- as.data.frame(matrix(data=NA,nrow=Duration,ncol=23))
  colnames(out) <- c("Date","StandAge","Year","StemNo","WF","WR","WS","avDBH","Height",
                     "StandVol","volWCer","volWCdr","NPP","NEE","LAI","Evapotranspiration","AvStemMass","BasArea","WSext",
                     "StandVol_loss", "VolProduction_tot", "DeepPercolation","RunOff")
  out[1,1] <- as.character(date)
  out[1,2:12] <- as.numeric(c(year,StandAge,StemNo,WF,WR,WS,avDBH,Height,StandVol,volWer,volWdr))
  
  #Starting Date of the climate data
  selectClimate <- as.numeric(which(date==climate[,1]))
  climate <- climate[selectClimate:length(climate[,1]),]
  
  ###########################################################
  #Do daily calculations
  ###########################################################
  for (day in 2:Duration){
    
    #next day of date
    date <- date+1
    
    #year
    year <- as.numeric(format(date,"%Y"))
    
    #Solar Radiation
    SolarRad <- climate$SolarRad[day]
    
    #Temperature data 
    Tav <- (climate$Tmax[day]+climate$Tmin[day])/2
    
    #VPD - mean day-time VPD (vapour pressure deficit)
    #gets daily "mean" VPD in mBar - based on daily max and min temperatures only
    #VPD is the difference (deficit) between the amount of moisture in the air and how much moisture the air can hold when it is saturated. 
    VPDx <- 6.1078 * exp(17.269 * climate$Tmax[day] / (237.3 + climate$Tmax[day]))
    VPDn <-  6.1078 * exp(17.269 * climate$Tmin[day] / (237.3 + climate$Tmin[day]))
    VPD <-  (VPDx + VPDn) / 2 #mean day-time VPD (vapour pressure deficit)
    
    #DayLength
    #gets fraction of day when sun is "up"
    dayofyr <- as.numeric(strftime(date,format="%j"))
    SLAt <- sin(round(pi,7) *lat / 180) #correct for rounding differences in pi
    cLat <- cos(round(pi,7) * lat / 180)
    sinDec <- 0.4 * sin(0.0172 * (dayofyr - 80))
    cosH0 <- -sinDec * SLAt / (cLat * sqrt (1- (sinDec)^2))
    if (cosH0 > 1) {
      getDayLength <- 0
    } else if (cosH0 < -1) {
      getDayLength <- 1
    } else {getDayLength <- acos(cosH0)/pi
    }
    
    #calculate seconds of the day when the sun is up
    DayLength <- 86400 * getDayLength
    
    #Precipitation
    Rain <- climate$Rain[day]
    
    #Determine various environmental modifiers
    
    #Temperature modifier (fT)
    if(Tav <= Tmin | Tav >= Tmax){
      fT <- 0
    }else{
      fT <- ((Tav - Tmin) / (Topt - Tmin)) * ((Tmax - Tav) / (Tmax - Topt)) ^ ((Tmax - Topt) / (Topt - Tmin))
    }
    
    #VPD modifier (fD)
    fVPD <- exp(-CoeffCond * VPD)
    
    #soil water modifier (fSW)
    MoistRatio <- min((erASW - erASWwp) / (erASWfc - erASWwp),1)
    fSW <- 1 / (1 + ((1 - MoistRatio) / SWconst) ^ SWpower)
    
    #Soil nutrition modifier (fN)
    if (fNn == 0) {
      fNutr <- 1 
    } else {
      fNutr <- 1 - (1 - fN0) * (1 - FR) ^ fNn
    }
    
    #Frost modifier (fF): if mean temperature < 0
    if(Tav<=0) {fFrost <- 1 - kF} else {fFrost <- 1} 
    
    #CO2 modifiers (fCaplpha/fCg)
    #CO2 added as functions for historical, RCP2.6, RCP8.5 (for now)
    if(CO2Concentration == "Historical") {CO2 <- CO2HistEq(year)}
    if(CO2Concentration == "RCP2.6") {CO2 <- CO2rcp26Eq(year)}
    if(CO2Concentration == "RCP4.5") {CO2 <- CO2rcp45Eq(year)}
    if(CO2Concentration == "RCP8.5") {CO2 <- CO2rcp85Eq(year)}
    if(class(CO2Concentration) == "numeric") {CO2 <- CO2Concentration}
    fCalpha <- fCalphax * CO2 / (350 * (fCalphax - 1) + CO2)
    fCg <- fCg0 / (1 + (fCg0 - 1) * CO2 / 350)
    
    #Age modifier (fage)
    if (nAge == 0) fAge <- 1
    RelAge <- StandAge / MaxAge
    fAge <- (1 / (1 + (RelAge / rAge) ^ nAge))
    
    #Physiological Modifier (PhysMod)
    PhysMod <- min(fVPD, fSW) * fAge 
    
    #Determine gross and net biomass production
    
    #Canopy cover and light interception (CanCover, lightIntcptn)
    if ((fullCanAge > 0) & (StandAge < fullCanAge)) {
      CanCover <- (StandAge + 0.01) / fullCanAge 
    } else {
      CanCover <- 1
    }
    lightIntcptn <- (1 - (exp(-k * LAI / CanCover)))
    
    #Calculate NPP
    alphaC <- alphaCx * fNutr * fT * fFrost * fCalpha * PhysMod  
    epsilon <- gDM_mol * molPAR_MJ * alphaC 
    RAD <- SolarRad
    RADint <- RAD * lightIntcptn #* CanCover
    RADsoil <- RAD - RADint
    GPP <- epsilon * RADint / 100  #convert into GPP/ha #tDM/ha day
    NPP <- GPP * y  
    
    #calculate canopy conductance and transpiration
    #CanCond <- getConductance(LAI,PhysMod,fCg)
    if (LAI <= LAIgcx) {
      gC <- MinCond + (MaxCond - MinCond) * LAI / LAIgcx
    } else {
      gC <- MaxCond
    }
    CanCond <-  gC * PhysMod * fCg
    if (CanCond == 0) CanCond <- 0.0001
    
    #Penman-Monteith equation for computing canopy transpiration
    #in kg/m2/day, which is converted to mm/day.
    netRad <- Qa + Qb * (RADint * 10 ^ 6 / DayLength)  #Q in MJ/m2/day --> W/m2
    defTerm <- rhoAir * lambda * (VPDconv * VPD) * BLcond
    div <-  CanCond * (1 + e20) + BLcond
    Etransp <- CanCond * (e20 * netRad + defTerm) / div      #in J/m2/s
    getTransp <- Etransp / lambda * DayLength    #converted to kg/m2/day
    Transp <- getTransp
    if(Transp < 0) Transp <- 0 #error control
    
    #rainfall interception
    if (LAImaxIntcptn > 0) {
      fracRainIntcptn <- MaxIntcptn * min(1, LAI / LAImaxIntcptn) 
    } else {
      fracRainIntcptn <- MaxIntcptn
    }
    
    #Snow Routine:
    if(Tav<=0 & Rain>0){
      Snow <- Rain
      Rain <- 0
      SnowInt <- Snow * fracRainIntcptn
      Snow <- Snow - SnowInt
      aSnow <- aSnow + Snow #Snow accumulation
      aSnowInt <- aSnowInt + SnowInt #interception snow accumulation
    }
    #Snowmelt:
    SnowMelt <- 0
    SnowIntMelt <- 0
    if(Tav>0 & aSnow > 0){
      SnowMelt <- min(SnowmMeltFactor*Tav,aSnow)
      SnowIntMelt <- min(SnowmMeltFactor*Tav,aSnowInt)
      aSnow <- aSnow - SnowMelt
      aSnowInt <- aSnowInt - SnowIntMelt
    }
    
    #Start snow period again without snow (1st of June)
    if(dayofyr==152){
      aSnow <- 0
      aSnowInt <- 0}
    
    #Rain Interception
    RainIntcptn <- Rain * fracRainIntcptn
    dayTF <- Rain - RainIntcptn #daily through fall
    
    #water balance: New sub-model
    #Set initials 0
    SupIrrig <- 0
    RunOff <- 0
    DP <- 0
    erPerc <- 0
    drPerc <- 0
    erUpflow <- 0
    drUpflow <- 0
    #Infiltration water
    pIW <- dayTF + pooledSW + SnowMelt #potential infiltration, later add snow melt
    if(pIW > maxInf){#infiltration cannot exceed max daily infiltration 
      IW <- maxInf
      RunOff <- pIW - IW
    } else {
      IW <- pIW
    }
    #Infiltration
    erASW <- erASW + IW
    excessSW <- max(erASW - erASWsat, 0)
    erASW <- erASW - excessSW
    pooledSW <- poolFractn * excessSW
    RunOff <- RunOff + (1 - poolFractn) * excessSW 
    #Transpiration
    EvapTransp <- min(erASW-erASWwp, Transp) #ET can not exceed erASW Wilting Point
    if(EvapTransp<0) EvapTransp <- 0 #error correction if erASW < erASWwp
    erASW <-  erASW - EvapTransp
    #Soil Evaporation
    if(erASW > erASWres){
      gSoil <- maxgSoil * ((erASW - erASWres) / (erASWsat - erASWres))  #0.0005 = maximum soil conductance, soil conductance
    } else{
      gSoil <- 0
    }
    #Penman-Monteith equation for computing soil evaporation
    #in kg/m2/day, which is converted to mm/day.
    netRad <- Qb * (RADsoil * 10 ^ 6 / DayLength)  #Q in MJ/m2/day --> W/m2, no Qa because all radiation on soil
    defTerm <- rhoAir * lambda * (VPDconv * VPD) * gAS
    div <-  gSoil * (1 + e20) + gAS
    ESoil <- gSoil * (e20 * netRad + defTerm) / div      #in J/m2/s
    getEvap <- ESoil / lambda * DayLength    #converted to kg/m2/day
    SoilEvap <- getEvap
    if(SoilEvap < 0 | Tav < 0 ) SoilEvap <- 0 #error control & if temperature <0°C no soil evaporation
    SoilEvap <- min(erASW-erASWres,SoilEvap) #Soil evaporation cannot exceed residual water content
    erASW <- erASW - SoilEvap
    
    #Add to Evapotranspiration
    EvapTransp <- EvapTransp + RainIntcptn + SnowIntMelt + SoilEvap
    
    #Percolation & Upflow calculations
    #unsaturated Percolation
    volWer <- erASW / tvER
    if(volWer <= volRes){
      WCer <- ((volRes * 1.0005) - volRes) / (volSat - volRes) #Water content cant exceed a fraction (0.05%) of residual water, needed for numerical solutions
    } else{
      WCer <- min((volWer - volRes) / (volSat - volRes), 1) #Water content ER
    }
    pMatricER <- (1 / VGalpha) * (((WCer ^ (-1 * VGn / (VGn - 1))) - 1) ^ (1 / VGn)) #matric potential ER, positive (after Van-Genuchten) but actually negative potential
    volWdr <- drASW / tvDR #Deep Root Zone
    if(volWdr <= volRes){
      WCdr <- ((volRes * 1.0005) - volRes) / ((drASWsat / tvDR) - volRes)
    } else{
      WCdr <- min((volWdr - volRes) / ((drASWsat / tvDR) - volRes), 1) 
    }
    pMatricDR <- (1 / VGalpha) * (((WCdr ^ (-1 * VGn / (VGn - 1))) - 1) ^ (1 / VGn))
    #From ER to DR percolation
    if(erASW >= erASWfc){
      dPm <- pMatricDR - pMatricER #delta matric potential > 0: faster Perc, dpM < 0: slower Perc
      if((dPm / depthER) < -1){
        dPm <- -1 * depthER #below -1: Upflow, because pMatric is higher than pGravity
      }
      unsatkF <- kS * (WCer ^ 0.5) * (1 - (1 - WCer ^ (VGn / (VGn - 1))) ^ (1 - 1 / VGn)) ^ 2 #relative hydraulic conductivity
      erPerc <- unsatkF * (dPm / depthER + 1) #Darcy-Buckingham-law; 2.part: water-potential gradient
      erPerc <- erPerc * 1000 #convert to mm/day
      erPerc <- min(erPerc, erASW - erASWfc)
      erASW <- erASW - erPerc
      drASW <- drASW + erPerc
    }
    #From DR to ER upflow
    if(erASW < erASWfc & drASW >= drASWfc & ((pMatricDR - pMatricER) / depthER) < -1){ #below -1 means upflow because pMatric higher than pGravity
      dPm <- pMatricDR - pMatricER
      unsatkF <- kS * (WCer ^ 0.5) * (1 - (1 - WCer ^ (VGn / (VGn - 1))) ^ (1 - 1 / VGn)) ^ 2
      erUpflow <- unsatkF * abs(dPm / depthER + 1)
      erUpflow <- erUpflow * 1000 
      erUpflow <- min(erUpflow, drASW - drASWfc) #Upflow cant exceed amount of SW under Field Capacity of DR-Zone
      erASW <- erASW + erUpflow
      drASW <- drASW - erUpflow
    }
    #DR to GW saturation flow
    if(drASW > drASWsat){
      drPerc <- drASW - drASWsat
      drASW <- drASW - drPerc
      DP <- drPerc #deep percolation ~ groundwater recharge
    }
    #Deep Root Percolation
    volWdr <- drASW / tvDR
    if(volWdr <= volRes){
      WCdr <- ((volRes * 1.0005) - volRes) / ((drASWsat / tvDR) - volRes)
    } else{
      WCdr <- min((volWdr - volRes) / ((drASWsat / tvDR) - volRes), 1)
    }
    pMatricDR <- (1 / VGalpha) * (((WCdr ^ (-1 * VGn / (VGn - 1))) - 1) ^ (1 / VGn))
    #DR to GW unsaturated flow
    if(drASW >= drASWfc & ((0 - pMatricDR) / depthDR) > -1){ #pMatricGW = 0, if <-1 then pMtric higher than pGravity
      dPm <- 0 - pMatricDR #delta matric potential DR
      unsatkF <- kSdr * (WCdr ^ 0.5) * (1 - (1 - WCdr ^ (VGn / (VGn - 1))) ^ (1 - 1 / VGn)) ^ 2
      drPerc <- unsatkF * (dPm / depthDR + 1) #Darcy-Buckingham-law
      drPerc <- drPerc * 1000 
      drPerc <- min(drPerc, drASW - drASWfc)
      drASW <- drASW - drPerc
      DP <- DP + drPerc
    }
    #GW to DR upflow
    if(drASW < drASWfc & ((0 - pMatricDR) / depthDR) < -1){
      dPm <- 0 - pMatricDR
      unsatkF <- kSdr * (WCdr ^ 0.5) * (1 - (1 - WCdr ^ (VGn / (VGn - 1))) ^ (1 - 1 / VGn)) ^ 2
      drUpflow <- unsatkF * Abs(dPm / depthDR + 1) #Abs beacuse negative (meaning Upflow) but for calculation needed positive
      drUpflow <- drUpflow * 1000 
      drASW <- drASW + drUpflow
      DP <- DP - drUpflow
    }
    #Error Control
    if(DP < 0) DP <- 0
    
    #Water Content Update
    volWer <- erASW/tvER
    volWdr <- drASW/tvDR
    
    #correct for actual ET
    if(Transp + RainIntcptn + SnowIntMelt + SoilEvap > 0){
      TranspScaleFactor <- EvapTransp / (Transp + RainIntcptn + SnowIntMelt + SoilEvap)
    }else{
      TranspScaleFactor <- 0
    }
    GPP <- TranspScaleFactor * GPP
    NPP <- TranspScaleFactor * NPP
    if(EvapTransp>0) WUE <- 100 * NPP / EvapTransp
    
    #Calculation of Net Ecosystem Production & Net Ecosystem Exchange (Meyer et al. 2018)
    #Temperature
    Rhsoil <- Rh10soil * Csoil * Q10 ^ (((0.0251 * Tav ^ 2 + 0.4079 * Tav + 1.0482) - 10) / 10)
    #  
    NEP <- NPP - Rhsoil #calculation of NEP (T DM/ha) 
    NEE <- -NEP                     
    
    ###############################################################################
    #Determine biomass increments and losses
    #calculate partitioning coefficients
    m <- m0 + (1 - m0) * FR
    pFS <- pfsConst * avDBH ^ pfsPower 
    pR <- pRx * pRn / (pRn + (pRx - pRn) * PhysMod * m)
    pS <- (1 - pR) / (1 + pFS) 
    pF <- 1 - pR - pS
    
    #calculate biomass increments
    incrWF <- NPP * pF
    incrWR <- NPP * pR
    incrWS <- NPP * pS
    #incrWF + incrWR + incrWS == NPP
    
    #calculate litterfall & root turnover -
    lossWF <- gammaF/30 * WF
    lossWR <- gammaR/30 * WR
    
    #Calculate end-of-month biomass
    WF <- WF + incrWF - lossWF
    WR <- WR + incrWR - lossWR
    WS <- WS + incrWS
    TotalW <- WF + WR + WS
    
    #Leaf fall
    MonthOneDayBefore <- currentMonth #needed for end of month calculations
    currentMonth <- as.numeric(format(as.Date(date,format="%d/%m/%Y"),"%m"))
    if(leaffall>0){
      currentDayMonth <- format(as.Date(date,format="%d/%m/%Y"),"%d-%m")
      if(currentDayMonth==paste0("01-",leaffall)) WFprior <- WF
      if(currentMonth==leaffall) WF <- max(WF - WFprior/31,0) #decrease dynamically over the month: end of leaffall no leaves
      if(currentMonth==leaffall+1) WF <- 0
      if(currentMonth==leafgrow) WF <- min(WF + WFprior/31,WFprior) #grow dynamically over the whole month: end of leafgrow WFprior 
    }
    
    #Update tree and stand data at the end of this time period,
    #taking mortality, thinning or defoliation into account
    
    #Age: daily age count, accounting for leap year
    lengthYear <- as.numeric(strftime(as.Date(paste0("31-12-",year),"%d-%m-%Y"),format = "%j"))
    StandAge <- StandAge + 1/lengthYear
    
    
    #Bark beetle attack, Thinning & Mortality: At the end of each month
    if(MonthOneDayBefore!=currentMonth){
      #Bark Beetle Attack, from Meyer et al. (2017) (could think to start it in april/may)
      if(StandAge >= attackAge & attackAge > 0){
        gammaNattack <- gammaN0attack * exp(-((StandAge - attackAge) * 12)/attackTime)
        delStems <- gammaNattack * StemNo / 100
        WF <- WF - mF * delStems * (WF / StemNo)
        WR <- WR - mR * delStems * (WR / StemNo)
        WS<- WS - mS * delStems * (WS / StemNo)
        StemNo <- StemNo - delStems
        mortality <- mortality + delStems
      }
      
      #Perform any thinning events: performed at the beginning of the year
      nThin <- as.numeric(length(thinAges))
      WSext <- 0
      if (thinEventNo <= nThin)  {
        if (StandAge >= thinAges[thinEventNo]) {
          if (StemNo > thinVals[thinEventNo]) { 
            delN <- (StemNo - thinVals[thinEventNo]) / StemNo
            StemNo <- StemNo * (1 - delN)
            WF <- WF * (1 - delN * thinWF[thinEventNo])
            WFprior <- WFprior * (1 - delN * thinWF[thinEventNo])
            WR <- WR * (1 - delN * thinWR[thinEventNo])
            WSext <- WS * delN * thinWS[thinEventNo]
            WS <- WS * (1 - delN * thinWS[thinEventNo])
          }
          thinEventNo <- thinEventNo + 1
        }
      }
      
      #Calculate age and stress-related mortality
      delStems <- 0
      mortality <- 0
      WSMort <- 0
      if (tgammaN != 0) {
        gammaN <- gammaN1 + (gammaN0 - gammaN1)*exp(-log(2)*(StandAge/tgammaN)^ngammaN)
      } else {
        gammaN <- gammaN1
      } 
      if (gammaN > 0) {
        delStems <- gammaN * StemNo / 12 / 100
        WF <- WF - mF * delStems * (WF / StemNo)
        WFprior <- WFprior - mF * delStems * (WFprior / StemNo)
        WR <- WR - mR * delStems * (WR / StemNo)
        WSmort <- mS * delStems * (WS / StemNo)
        WS <- WS - mS * delStems * (WS / StemNo)
        StemNo <- StemNo - delStems
        mortality <- mortality + delStems
      }
      
      #Calculate self-thinning mortality
      selfThin <- 0
      WSselfThin <- 0
      wSmax <- wSx1000 * (1000 / StemNo) ^ thinPower #(1000 / StemNo) ^ thinPower can be seen as a factor/ empirical relationship determining how soon selfthinnning starts
      AvStemMass <- WS * 1000 / StemNo #transform Ws into kg/tree
      delStems <- 0
      if (wSmax < AvStemMass) {
        accuracy <- 1 / 1000
        n <- StemNo / 1000
        x1 <- 1000 * mS * WS / StemNo
        for (i in 1:5){
          x2 <- wSx1000 * n ^ (1-thinPower)  
          fN <- x2 - x1 * n - (1 - mS) * WS
          dfN <- (1 - thinPower) * x2 / n - x1
          dN <- -fN / dfN
          n <- n + dN
          if (abs(dN)<=accuracy) break
        }
        delStems <- StemNo - 1000 * n
        
        WF <- WF - mF * delStems * (WF / StemNo)
        WFprior <- WFprior - mF * delStems * (WFprior / StemNo)
        WR <- WR - mR * delStems * (WR / StemNo)
        WSselfThin <- mS * delStems * (WS / StemNo)
        WS <- WS - mS * delStems * (WS / StemNo)
        StemNo <- StemNo - delStems
        wSmax <- wSx1000 * (1000 / StemNo) ^ thinPower
        AvStemMass <- WS * 1000 / StemNo
        selfThin <- selfThin + delStems
      }
    }
    
    ##################################################################
    #update age dependent factors
    #SLA = expF(StandAge, SLA0, SLA1, tSLA, 2)
    if (tSLA != 0) {
      SLA <- SLA1 + (SLA0 - SLA1)*exp(-log(2)*(StandAge/tSLA)^2) 
    } else {
      SLA <- SLA1
    }
    
    #branch and bark fraction
    if (tBB != 0) {
      fracBB <- fracBB1 + (fracBB0 - fracBB1)*exp(-log(2)*(StandAge/tBB)^1) 
    } else {
      fracBB <- fracBB1
    }
    
    #Density = expF(StandAge, rho0, rho1, tRho, 1)
    if (tRho != 0) {
      Density <- rho1 + (rho0 - rho1)*exp(-log(2)*(StandAge/tRho)^1) 
    } else {
      Density <- rho1
    }
    
    #gammaF <- gammaFoliage(StandAge)
    if (tgammaF * gammaF1 == 0) {
      gammaF <- gammaF1 
    } else {
      kgammaF <- 12 * log(1 + gammaF1 / gammaF0) / tgammaF
      gammaF <- gammaF1 * gammaF0 / (gammaF0 + (gammaF1 - gammaF0) * exp(-kgammaF * StandAge))
    }
    
    #update stand characteristics
    LAI <- WF * SLA * 0.1
    avDBH <- (AvStemMass / aWs) ^ (1 / nWs) 
    BasArea <- ((avDBH / 200) ^ 2) * pi * StemNo
    if(HeightEquation==1) Height <-  aH * avDBH ^ nHB * StemNo ^ nHC
    if(HeightEquation==2) Height <- 1.3 + aH * exp(-nHB/avDBH) + nHC * Density * avDBH #Michajlow-Schumacher (Forrester et. al, 2021)
    StandVol_loss <- 0
    if(SVEquation==1) StandVol <- aV * (avDBH ^ nVB) * (Height ^ nVH) * StemNo #equation with parameters after Forrester et al. 2021
    if(SVEquation==2) StandVol <- WS * (1 - fracBB) / Density #3PG original equation
    if(SVEquation==3) StandVol <- 0.7 * ((avDBH/200)^2*pi) * Height * StemNo #simple volume equation
    if(StandVol < oldV) StandVol_loss <- oldV - StandVol
    oldV <- StandVol
    
    #accumulated volume: Production
    accV <- accV+StandVol_loss
    VolProduction_tot <- StandVol + accV
    
    #MAI
    #MAI <- (StandVol + StandVol_ext)/StandAge
    
    #WF in kg per tree
    #WFtree <- WF *1000 /StemNo
    out[day,1] <- as.character(date)
    out[day,2:23] <- as.numeric(c(year,StandAge,StemNo,WF,WR,WS,avDBH,Height,StandVol,volWer,volWdr,NPP,NEE,LAI,EvapTransp,
                                  AvStemMass,BasArea,WSext,StandVol_loss, VolProduction_tot,DP,RunOff))
  }
  
  ###########################################################
  #Create Final output data
  ###########################################################
  out$Date <- as.Date(out$Date,format="%Y-%m-%d")
  return(out)
}
