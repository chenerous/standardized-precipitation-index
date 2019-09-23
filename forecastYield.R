####  forecastYield.R 
####
####  test models for yields by their forecasting ability
####  build the models using 25 years of data
####  then test on the last 5 years of data
####
####  compare a stepwise linear regression model with a random forest model 
####
#### inputs:  
####  output of the SPI code
####  a US county shapefile 
####  NASS rainfed corn yieldsfrom https://quickstats.nass.usda.gov
####
#### shapefile from:
#### https://catalog.data.gov/dataset/tiger-line-shapefile-2016-nation-u-s-current-county-and-equivalent-national-shapefile/resource/fa774c9d-a098-4792-bfd4-94c7caa190b6
####
#### libraries required: sf, MASS, randomForest
####
#### 

library(sf)
library(MASS)
library(randomForest)


###########################################################################
#### data reading and organization
#### at the end, we have a set of regressors in zplot_ordered (spi by month by county)
#### and outcomes in yields_ordered (county corn yield)
#### but yields have not yet been detrended
###########################################################################

#### read SPI data
load("spi_shiny_app/data/spi/SPI_out_3month.RData")

### pick out which estimate you want to use z, z_emp, z_fit
zplot <- spi_results$z
Lon = spi_results$Lon
Lat = spi_results$Lat
monthList = spi_results$monthList
yearList = spi_results$yearList
yearSeq = sort(unique(yearList))
numYears = length(yearSeq)


rm(spi_results)
for (i in 1:30) gc()

#### read below for bigger shapefile with lat lon data
counties <- st_read("spi_shiny_app/data/us_counties_large/tl_2016_us_county.shp",stringsAsFactors=FALSE)

#### assign county SPI value to closest cell of z 
#### (calculating a county average might be nice to add later, but averaging effect might actually be detrimental to model)

#### create map of counties to cells
#### reorder zplot data
#### might be faster way to do this (to avoid loop), but only 3108 counties
celltocountymap <- array(0,c(nrow(counties),3))
zplot_ordered <- array(0,c(nrow(counties),dim(zplot)[3]))
for (icounty in 1:nrow(counties)){
  i=which(abs(as.numeric(counties$INTPTLON[icounty])-Lon) == min(abs(as.numeric(counties$INTPTLON[icounty])-Lon)))
  j=which(abs(as.numeric(counties$INTPTLAT[icounty])-Lat) == min(abs(as.numeric(counties$INTPTLAT[icounty])-Lat)))
  celltocountymap[icounty,1]=icounty
  celltocountymap[icounty,2]=i
  celltocountymap[icounty,3]=j
  zplot_ordered[icounty,]=zplot[i,j,]
}
#rm(zplot)
#for (i in 1:30) gc()

stateList = counties$STATEFP
countyList = counties$COUNTYFP
IDList = counties$GEOID
countyLatList = counties$INTPTLAT
countyLonList = counties$INTPTLON

rm(counties)
for (i in 1:30) gc()

yields <- read.csv("../NASS_corn_1989_2008.csv")
yields <- rbind(yields,read.csv("../NASS_corn_2009_2018.csv"))
yields_ordered = array(NA,c(dim(zplot_ordered)[1],numYears))
for (i in 1:length(stateList)){
  i0 <-which(yields$State.ANSI == as.numeric(stateList[i]) & yields$County.ANSI == as.numeric(countyList[i]))
  j0 <- match(yields$Year[i0],yearSeq)
  yields_ordered[i,j0]=yields$Value[i0]
}

#######################################################
#### detrend yields (linear trends only right now!)
######################################################
for (i in 1:length(countyList)){
  if (!all(is.na(yields_ordered[i,]))){
        
    model <- lm(yields_ordered[i,]~seq(1,numYears),na.action=na.exclude)
    trend <- predict(model)
    yields_ordered[i,]=yields_ordered[i,]-trend
  }
}



#######################################################
#### model A: linear stepwise regression
#### only use weather data from May to October (generic crop season)
#### this means 6 regressors per yield value (25 years of yields)
######################################################

modYears=25
predictYears =5
predictedYields <- array(NA,c(dim(zplot_ordered)[1],numYears))

#### for each county, build a model, then test it
for (i in 1:length(countyList)){
  if (length(which(!is.na(yields_ordered[i,]) & abs(yields_ordered[i,]) > .1)) > 1){
    #### create model data frame
    y=yields_ordered[i,1:modYears]
    x1=zplot_ordered[i,which(monthList==5)[1:modYears]]
    x2=zplot_ordered[i,which(monthList==6)[1:modYears]]
    x3=zplot_ordered[i,which(monthList==7)[1:modYears]]
    x4=zplot_ordered[i,which(monthList==8)[1:modYears]]
    x5=zplot_ordered[i,which(monthList==9)[1:modYears]]
    x6=zplot_ordered[i,which(monthList==10)[1:modYears]]
    mydata=data.frame(y,x1,x2,x3,x4,x5,x6)
    #### use stepwise linear regression to pick out important regressors
    model <- lm(y~x1+x2+x3+x4+x5+x6,data=mydata,na.action=na.exclude)
    result <- try(step<- stepAIC(model,direction='both'))
    if (class(result) != "try-error"){
      formula <- formula(step)}
    #### create pred data frame
    #y=yields_ordered[i,(numYears-predictYears+1):numYears]
    y=yields_ordered[i,]
    if (!all(is.na(y))){
      # x1=zplot_ordered[i,which(monthList==5)[(numYears-predictYears+1):numYears]]
      # x2=zplot_ordered[i,which(monthList==6)[(numYears-predictYears+1):numYears]]
      # x3=zplot_ordered[i,which(monthList==7)[(numYears-predictYears+1):numYears]]
      # x4=zplot_ordered[i,which(monthList==8)[(numYears-predictYears+1):numYears]]
      # x5=zplot_ordered[i,which(monthList==9)[(numYears-predictYears+1):numYears]]
      # x6=zplot_ordered[i,which(monthList==10)[(numYears-predictYears+1):numYears]]
      x1=zplot_ordered[i,which(monthList==5)]
      x2=zplot_ordered[i,which(monthList==6)]
      x3=zplot_ordered[i,which(monthList==7)]
      x4=zplot_ordered[i,which(monthList==8)]
      x5=zplot_ordered[i,which(monthList==9)]
      x6=zplot_ordered[i,which(monthList==10)]
      #### pull out only important regressors and predict
      if (class(result) != "try-error"){
        if (dim(step$model)[2] > 1){
          varnames <- names(step$model)[2:length(names(step$model))]
          newdata <- data.frame(y)
          for (j in 1:length(varnames)){
            newdata <- cbind(newdata,get(varnames[j]))
          }
          names(newdata)[2:(1+length(varnames))]=varnames
          model <- lm(formula,newdata,na.action=na.exclude)
          newdata <- data.frame(y)
          for (j in 1:length(varnames)){
            newdata <- cbind(newdata,get(varnames[j]))
          }
          newdata <- newdata[2:(1+length(varnames))]
          names(newdata)=varnames
        }else{
          newdata=data.frame(y,x1,x2,x3,x4,x5,x6)
        }
      }
      predictedYields[i,]=predict(model,newdata)
    }
  }
}

#######################################################
#### model B: random forest
######################################################

######################################################
### below is a random forest using only one county's data
######################################################
# #### for each county, build a model, then test it
# #### ordered data is 7 features of 30 observations
# for (i in 1:length(countyList)){
#   if (!all(is.na(yields_ordered[i,]))){
#     #### order data
#     yield=yields_ordered[i,]
#     month5=zplot_ordered[i,which(monthList==5)]
#     month6=zplot_ordered[i,which(monthList==6)]
#     month7=zplot_ordered[i,which(monthList==7)]
#     month8=zplot_ordered[i,which(monthList==8)]
#     month9=zplot_ordered[i,which(monthList==9)]
#     month10=zplot_ordered[i,which(monthList==10)]
#     County.Corn = data.frame(yield,month5,month6,month7,month8,month9,month10)  
#     #### fix below (need to deal with NAs)
#     set.seed(1)
#     bag.County.Corn=randomForest(yield~.,data=County.Corn,subset=seq(1,25),mtry=6,importance=TRUE)
#     bag.County.Corn
#     yhat.bag = predict(bag.County.Corn,newdata=County.Corn)
#     plot(yhat.bag, County.Corn$yield)
#     points(yhat.bag[25:30], County.Corn$yield[25:30],col='red')
#     abline(0,1)
#   }
# }

######################################################
### below is a random forest using only all county data
######################################################
#### ordered data features:  yields, 6 1-month spi values, county, year  
yield=c()
month5=c()
month6=c()
month7=c()
month8=c()
month9=c()
month10=c()
county=c()
year=c()
for (i in 1:length(countyList)){
  if (!all(is.na(yields_ordered[i,]))){
    for (j in 1:length(yearSeq)){
      if (!is.na(yields_ordered[i,j])){
        #### order data
        yield=c(yield,yields_ordered[i,j])
        month5=c(month5,zplot_ordered[i,which(monthList==5 & yearList == yearSeq[j])])
        month6=c(month6,zplot_ordered[i,which(monthList==6 & yearList == yearSeq[j])])
        month7=c(month7,zplot_ordered[i,which(monthList==7 & yearList == yearSeq[j])])
        month8=c(month8,zplot_ordered[i,which(monthList==8 & yearList == yearSeq[j])])
        month9=c(month9,zplot_ordered[i,which(monthList==9 & yearList == yearSeq[j])])
        month10=c(month10,zplot_ordered[i,which(monthList==10 & yearList == yearSeq[j])])
        county=c(county,as.numeric(IDList[i]))
        year=c(year,yearSeq[j])
      }
    }
  }
}
County.Corn.unordered = data.frame(yield,month5,month6,month7,month8,month9,month10,county,year)  
i0 <- order(year)
County.Corn <- County.Corn.unordered[i0,]

### do random forest on all data before 2013
print('start random forest')
i0 <- which(County.Corn$year <= 2013)
set.seed(1)
bag.County.Corn=randomForest(yield~.,data=County.Corn,subset=i0,mtry=6,importance=TRUE)
print('end random forest')

### predict yields
yhat.bag = predict(bag.County.Corn,newdata=County.Corn)

#######################################################
#### compare results:  Mclean County, IL
###
### Mclean County is one of the largest corn producing counties in the US
### FIPS = 17113
######################################################
countyfips=17113
i0 <- which(as.numeric(IDList) == countyfips)
plot(yearSeq,yields_ordered[i0,],pch=16,cex=1,col='red',type='b',xlab='Year',ylab='Yield Relative to Normal (bu/ac)',main='Example of Yield Prediction Performance')
points(yearSeq,predictedYields[i0,],pch=16,type='b')

j0 <- which(County.Corn$county == countyfips)
points(yearSeq,yhat.bag[j0],col='blue',pch=16,type='b')
abline(v=2014)
text(2000,-60,c("Red points are observed yields"),cex=0.7)
text(2000,-65,c("Black points are regression prediction"),cex=0.7)
text(2000,-70,c("Blue points are random forest prediction"),cex=0.7)
text(2000,-75,c("Model is trained on data before 2014"),cex=0.7)
text(2000,-80,c("Model is predicting data after 2013"),cex=0.7)

#######################################################
#### compare results:  scatter
######################################################

par(mfrow=c(1,2))
i0 <- which(yearSeq >=2014)
plot(yields_ordered[,i0],predictedYields[,i0],pch=16,cex=0.5,xlim=range(-100,130),ylim=range(-100,130),xlab='Observed Yield Dev',ylab='Modeled Yield Dev',main='Stepwide Regeression')
abline(0,1,col='red')
j0 <- which(County.Corn$year >= 2014)
plot(County.Corn$yield[j0],yhat.bag[j0],col='blue',pch=16,cex=0.5,xlim=range(-100,130),ylim=range(-100,130),xlab='Observed Yield Dev',ylab='Modeled Yield Dev',main='Random Forest')
abline(0,1,col='red')

#######################################################
#### compare results:  explaining the results
######################################################

#bag.County.Corn
#varImpPlot(bag.County.Corn)
