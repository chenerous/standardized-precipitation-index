#####################################################################
#### title:  calcSPI.R
#### 
#### purpose:  calculate standardized precipitation index
####
#### add references
####
#### Input:  daily precipitation of 0.25 lat-lon grid in units of 0.1mm
#### Steps:  
#### 1. read in precip
#### 2. aggregate to month 
#### 3. calculate total precipitation for time scale of interest at monthly intervals
####    for example, if the time scale is 3 months, then the data should look like:
####        month 1 = NA
####        month 2 = NA
####        month 3 = total rain in months 1-3
####        month 4 = total rain in months 2-4
#####################################################################

#### read and eval config file
configFile <- "spi.config"
source("readConfig.R")


#### set up some stuff (values, arrays)
numYears = length(seq(yearFrom,yearTo))
numMonths=numYears*12
monthList<- rep(seq(1,12),length(numYears))
yearList<- c()
for (year in yearFrom:yearTo){yearList <- c(yearList,rep(year,12))}
nLon = 1440 # the gridding of the NOAA CPC gridded CONUS precip data 
nLat = 720 # the gridding of the NOAA CPC gridded CONUS precip data 
precData <- array(0,c(nLon,nLat,numMonths))
intervalData <- array(0,c(nLon,nLat,numMonths))
intervalRankData <- array(0,c(nLon,nLat,numMonths))

###########################################
# read precip data and aggregate to month
# precip data is stored in PrecData
###########################################

# get list of files
filelist <- list.files(dirPrecip)

# for each file 
for (filename0 in filelist){
  
  # pull out date from filename
  yearStamp = as.numeric(substr(filename0,37,40))
  monthStamp = as.numeric(substr(filename0,41,42))
  
  # create filename
  filename = paste(dirPrecip,filename0,sep="")
  
  # consider as gz or not (I think the CONUS data may all be non gz)
  if (length(grep(".gz",filename)) > 0) {
    con <- gzfile(filename,open="rb")
  } else {
    con <- file(filename,open="rb")
  }
  # read file
  precCONUS <- readBin(con,"numeric",size=4,endia="little",n=nLon*nLat)
  
  # put into array (not totally sure about order of lon and lat -- will need to check)
  precCONUS <- array(precCONUS,c(nLon,nLat)) #array(round(precCONUS),c(nLon,nLat))
  
  # this would be a good point to do some data checks or cleaning
  
  # close file
  close(con)
  
  # add to month data 
  i0 <- which(monthList == monthStamp & yearList == yearStamp)
  precData[,,i0] <- precData[,,i0]+precCONUS
}


###########################################
# calculate total precip for time interval
# put in intervalData
# rank the years for precip for each cell and month of year
###########################################

for (i in 1:(spi_interval-1)){
 intervalData[,,i] = NA 
}
for (i in spi_interval:numMonths){
  for (j in 1:spi_interval){
  intervalData[,,i]=intervalData[,,i]+precData[,,(i-(j-1))]
  }
}

for (iMonth in 1:12) {
  monthSequence=seq(iMonth,numMonths,by=12)
  intervalRankData[,,MonSeq]<- apply(intervalData[,,monthSequence],c(1,2),rank,na.last=FALSE)
}


###########################################
# calculate sum over years 
# calculate the probability of zero precip (q)
# estimate parameters in gamma function (ahat,bhat)
###########################################

sumZero <- array(0,c(nLon,nLat,12))
sumNonZero <- array(0,c(nLon,nLat,12))
sumLogPrec <- array(0,c(nLon,nLat,12))
ahat <- array(0,c(nLon,nLat,12))
bhat <- array(0,c(nLon,nLat,12))  

for (iMonth in 1:12) {
  monthSequence=seq(iMonth,numMonths,by=12)
  testIntervalData=intervalData[,,monthSequence]
  logIntPrec = log(testIntervalData)   
  sumPrec=apply(testIntervalData,c(1,2),sum,na.rm=TRUE)
  
  indZero <- testIntervalData == 0
  indPos <- testIntervalData > 0
  for (i in 1:numCells){
    sumZero[,,iMonth] = apply(indZero,c(1,2),sum,na.rm=TRUE)
    sumNonZero[,,iMonth] = apply(indPos,c(1,2),sum,na.rm=TRUE)
    #  I think I just need to sum the log of non-zero precip
    sumLogPrec[,,iMonth]=apply(logIntPrec[indPos],c(1,2),sum,na.rm=TRUE)
    #
  }
  avgSumPrec=sumPrec/sumNonZero[,,iMonth]
  
  bigA = log(avgSumPrec) - sumLogPrec[,,iMonth]/sumNonZero[,,iMonth] 
  ahat[,,iMonth]=1./(4.*bigA)*(1.+sqrt(1.+4.*bigA/3.))
  bhat[,,iMonth]=avgSumPrec/ahat[,,iMonth]
}

ProbZero=sumZero/(sumZero+sumNonZero)


###########################################
# do fit of gamma function with ahat and bhat and starting points
###########################################

afit <- ahat
bfit <- bhat
for (i in 1:nLon){
  for (j in 1:nLat){
    for (iMonth in 1:12) {
      monthSequence=seq(iMonth,numMonths,by=12)
      if (monthSequence[1] < spi_interval) {monthSequence <- monthSequence[2:length(monthSequence)]}
      IntPrec_fit=intervalData[i,j,monthSequence]
      ### below line needs checking
      IntPrec_fit <- IntPrec_fit[!is.na(IntPrec_fit)] 
      if (length(IntPrec_fit) >= 30) {
        # the code may not like doing the fit and some playing wih 'upper','lower', or other fitdistr inputs may be needed
        testFitParam <- fitdistr(IntPrec_fit,"gamma",list(shape=ahat[i,j,iMonth],scale=bhat[i,j,iMonth]))
        afit[i,j,iMonth] = testFitParam$estimate[1]
        bfit[i,j,iMonth] = testFitParam$estimate[2]
      }
    }
  }
}

###########################################
# calculate cumulative probability
# first incomplete gamma function (G(x))
# then tot prob with zeros (H(x))
###########################################

H <- array(0, c(nLon,nLat,numMonths))
tX <- array(0,c(maxMonths))
H_fit <- array(0, c(nLon,nLat,numMonths))
tX_fit <- array(0,c(numMonths))

for (i in 1:nLon){
  for (j in 1:nLat){
  for (iMonth in 1:12) {
    monthSequence=seq(iMonth,numMonths,by=12)
    tX[monthSequence]=IntPrec[i,j,monthSequence]/bhat[i,j,iMonth]
    tX_fit[monthSequence]=IntPrec[i,j,monthSequence]/bfit[i,j,iMonth]
    
    
    H[i,j,monthSequence]=pgamma(tX[monthSequence],ahat[i,iMonth])
    H[i,j,monthSequence] =ProbZero[i,j,iMonth]+(1.-ProbZero[i,j,iMonth])*H[i,j,monthSequence]
    H_fit[i,j,monthSequence]=pgamma(tX_fit[monthSequence],afit[i,j,iMonth])
    H_fit[i,j,monthSequence] =ProbZero[i,j,iMonth]+(1.-ProbZero[i,j,iMonth])*H_fit[i,j,monthSequence]
  }
}}

###########################################
# Do equiprobability transformation from H(x) to SPI
# via Abramowitz & Stegun 
# result is SPI in gridded array for every month and every year
###########################################

c0=2.515517
c1=0.802853
c2=0.010328
d1=1.432788
d2=0.189269
d3=0.001308

z <- array(0,c(nLon,nLat,numMonths))
z_emp <- array(0,c(nLon,nLat,numMonths))
tH <- array(0,c(nLon,nLat,numMonths))
z_fit <- array(0,c(nLon,nLat,numMonths))
tH_fit <- array(0,c(nLon,nLat,numMonths))

H[is.na(H)] <- 0
ind0 <- H <= .5 
ind1 <- (H > 0.5 & H <= 1.)
tH[ind0] = sqrt(log(1/H[ind0]^2))
tH[ind1] = sqrt(log(1/(1.-H[ind1])^2))

z=(tH-(c0+c1*tH+c2*tH^2)/(1+d1*tH+d2*tH^2+d3*tH^3))
z[ind0]=-z[ind0]

H_fit[is.na(H_fit)] <- 0
ind0 <- H_fit <= .5 
ind1 <- (H_fit > 0.5 & H_fit <= 1.)
tH_fit[ind0] = sqrt(log(1/H_fit[ind0]^2))
tH_fit[ind1] = sqrt(log(1/(1.-H_fit[ind1])^2))

z_fit=(tH_fit-(c0+c1*tH_fit+c2*tH_fit^2)/(1+d1*tH_fit+d2*tH_fit^2+d3*tH_fit^3))
z_fit[ind0]=-z_fit[ind0]

# calculated empirical values
for (i in 1:nLon){
  for (j in 1:nLat){
  for (iMonth in 1:12) {
    monthSequence=seq(iMonth,numMonths,by=12)
    # all the right SPI values, not in the right order
    
    z_emp[i,j,monthSequence]=seq(1,length(monthSequence))/length(monthSequence)
    if (monthSequence[1] < spi_interval) {     
      z_emp[i,j,monthSequence[2:length(monthSequence)]]=seq(1,length(monthSequence)-1)/(length(monthSequence)-1)
      z_emp[i,j,monthSequence[1]]=NA
    }
    
    z_emp[i,j,monthSequence]=qnorm(z_emp[i,j,monthSequence])
    # fix order somehow
    z_emp[i,j,monthSequence]=z_emp[i,j,iRankIntPrec[i,j,monthSequence]]
  }
} }

# cap everything to be between -3 and 3
z[z>3] <- 3
z[z< -3] <- -3
z_fit[z_fit>3] <- 3
z_fit[z_fit< -3] <- -3
z_emp[z_emp>3] <- 3
z_emp[z_emp< -3] <- -3


###########################################
# save results
###########################################


spi_results = list(z,z_emp,z_fit)
save(spi_results,file=fileSPIout)
