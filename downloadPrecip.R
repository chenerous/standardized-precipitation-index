#### title:  downloadPrecip.R
#### 
#### purpose:  download NOAA CPC gridded CONUS data from 
####  ftp://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_CONUS/
####
#### notes:
#### this code requires the RCurl library
#### to calculate the spi, need to have at least 30 years of data

#### read and eval config file
configFile <- "spi.config"
source("readConfig.R")

#### download sufficient data
yearList = seq(yearFrom,yearTo)
for (year in yearList){
  if (year <= 2006){
    url0 = "ftp://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_CONUS/V1.0/"
  }else{
    url0 = "ftp://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_CONUS/RT/"
  }
url=paste(url0,year,"/",sep="")
filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames = unlist(filenames)
for (filename in filenames) {
  download.file(paste(url, filename, sep = ""), paste(dirPrecip,filename,sep = ""))
}
}


