####  plotSPI_leaflet.R
####
####  create map of SPI using leaflet 
####
#### inputs:  
####  output of the SPI code
####  a US county shapefile 
####
#### shapefiles from:
#### https://catalog.data.gov/dataset/tiger-line-shapefile-2016-nation-u-s-current-county-and-equivalent-national-shapefile/resource/fa774c9d-a098-4792-bfd4-94c7caa190b6
#### https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
####
#### libraries required: leaflet, sf
####
#### 

library(leaflet)
library(sf)

#### read SPI data
load("SPI_out_3month.RData")

### pick out which estimate you want to use z, z_emp, z_fit
zplot <- spi_results$z
Lon = spi_results$Lon
Lat = spi_results$Lat

### pick out which month and year you want
month = 6
year = 2006

iMonth=which(spi_results$monthList == month & spi_results$yearList == year)

### read shapefile, project to WGS84

counties <- st_read("../cb_2018_us_county_20m/cb_2018_us_county_20m.shp",stringsAsFactors=FALSE)
counties <- st_transform(counties,crs ="+proj=longlat +datum=WGS84")

#### read below for bigger shapefile with lat lon data
counties_wlatlon <- st_read("../tl_2016_us_county/tl_2016_us_county.shp",stringsAsFactors=FALSE)

#### pull out only contiguous US
contig_counties <- subset(counties,!((counties$STATEFP) %in% c("02","03","07","14","15","43","52","60","66","69","72","78")))

#### match lat-lon
i0 <- match(contig_counties$GEOID,counties_wlatlon$GEOID)
contig_counties$Lat = counties_wlatlon$INTPTLAT[i0]
contig_counties$Lon = counties_wlatlon$INTPTLON[i0]

#### color counties by closest cell value of z
#### calculating a county average would be nice to add later
#### might be faster way to do this (to avoid loop), but only 3108 counties
contig_counties$SPI= 0
for (icounty in 1:nrow(contig_counties)){
  i=which(abs(as.numeric(contig_counties$Lon[icounty])-Lon) == min(abs(as.numeric(contig_counties$Lon[icounty])-Lon)))
  j=which(abs(as.numeric(contig_counties$Lat[icounty])-Lat) == min(abs(as.numeric(contig_counties$Lat[icounty])-Lat)))
  contig_counties$SPI[icounty]=plot[i,j,iMonth]
}

#### create leaflet
bins <- c(-3,-2,-1.5,-1,1,1.5,2,3)
pal <- colorBin("RdYlBu", domain = contig_counties$SPI, bins = bins)

m<-leaflet(contig_counties) %>%
  addTiles() %>%
  setView(-95, 38, zoom = 4) %>%
  addPolygons(label = ~paste(NAME,": ",round(SPI,2),sep=""),
              fillColor = ~pal(SPI),
              stroke=FALSE,
              smoothFactor = 0,
              fillOpacity = 0.5,
              highlightOptions = highlightOptions(color = "white", weight = 2,
                                                  bringToFront = TRUE)) %>%
  addLegend(pal = pal,
          values = ~SPI,
          opacity = 0.7,
          title =  paste(as.character(month),"/",as.character(year),' 3 month SPI',sep=""),
          position = "bottomright") 

#### above map has stroke=FALSE  => no lines between counties
#### to put those lines back in, delete stroke=FALSE and add
####  color = "#444444", weight = 1, 

