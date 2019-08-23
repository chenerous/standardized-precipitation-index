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
precData = array(0,c(nLon,nLat,numMonths))

#### read precip data and aggregate to month

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

