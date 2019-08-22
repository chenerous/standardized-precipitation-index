#### title:  downloadPrecip.R
#### 
#### purpose:  download NOAA CPC gridded CONUS data from 
####  ftp://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_CONUS/
####
#### this code requires the RCurl library

#### read and eval config file
configFile <- "spi.config"
source("readConfig.R")

#### edit to download sufficient data
url = "ftp://ftp.cpc.ncep.noaa.gov/precip/CPC_UNI_PRCP/GAUGE_CONUS/RT/2018/"
filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
filenames <- strsplit(filenames, "\n")
filenames = unlist(filenames)
for (filename in filenames[1]) {
  download.file(paste(url, filename, sep = ""), paste("CPC_precip_data/",filename,sep = ""))
}

