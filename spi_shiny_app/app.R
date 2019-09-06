####  Shiny app implementation of plotSPI_leaflet.R
####
####  create map of SPI using leaflet 
####
#### inputs:  
####  output of the SPI code
####  a US county shapefile 
####
#### shapefile from:
#### https://catalog.data.gov/dataset/tiger-line-shapefile-2016-nation-u-s-current-county-and-equivalent-national-shapefile/resource/fa774c9d-a098-4792-bfd4-94c7caa190b6
#### https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
####
#### libraries required: leaflet, sf
####
#### 

library(shiny)
library(leaflet)
library(sf)

###########################################################################
#### all the reading of data and processing that only needs to be done once
###########################################################################

#### read SPI data
load("data/spi/SPI_out_3month.RData")

### pick out which estimate you want to use z, z_emp, z_fit
zplot <- spi_results$z
Lon = spi_results$Lon
Lat = spi_results$Lat
monthList = spi_results$monthList
yearList = spi_results$yearList

rm(spi_results)
for (i in 1:30) gc()

### read shapefile, project to WGS84
counties <- st_read("data/us_counties_small/cb_2018_us_county_20m.shp",stringsAsFactors=FALSE)
counties <- st_transform(counties,crs ="+proj=longlat +datum=WGS84")

#### pull out only contiguous US
contig_counties <- subset(counties,!((counties$STATEFP) %in% c("02","03","07","14","15","43","52","60","66","69","72","78")))

rm(counties)
for (i in 1:30) gc()

#### read below for bigger shapefile with lat lon data
counties_wlatlon <- st_read("data/us_counties_large/tl_2016_us_county.shp",stringsAsFactors=FALSE)

#### match lat-lon
i0 <- match(contig_counties$GEOID,counties_wlatlon$GEOID)
contig_counties$Lat = counties_wlatlon$INTPTLAT[i0]
contig_counties$Lon = counties_wlatlon$INTPTLON[i0]

rm(counties_wlatlon)
for (i in 1:30) gc()

#### create dataset of county values for every month in spi time series, then only pick correct data set when rendering (below)

#### assign county SPI value to closest cell of z 
#### (calculating a county average might be nice to add later, but not sure averaging effect is worthwhile)

#### create map of counties to cells
#### reorder zplot data
#### might be faster way to do this (to avoid loop), but only 3108 counties
celltocountymap <- array(0,c(nrow(contig_counties),3))
zplot_ordered <- array(0,c(nrow(contig_counties),dim(zplot)[3]))
for (icounty in 1:nrow(contig_counties)){
  i=which(abs(as.numeric(contig_counties$Lon[icounty])-Lon) == min(abs(as.numeric(contig_counties$Lon[icounty])-Lon)))
  j=which(abs(as.numeric(contig_counties$Lat[icounty])-Lat) == min(abs(as.numeric(contig_counties$Lat[icounty])-Lat)))
  celltocountymap[icounty,1]=icounty
  celltocountymap[icounty,2]=i
  celltocountymap[icounty,3]=j
  zplot_ordered[icounty,]=zplot[i,j,]
}
rm(zplot)
for (i in 1:30) gc()

#### placeholder for SPI values in the contig_counties data set
contig_counties$SPI= 0

###########################################################################
#### ui
###########################################################################

ui <- fluidPage(
  
  titlePanel("3 month SPI"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("3 month Standardized Precipitation for counties in the contiguous United States, 1989-2018."),
      
      sliderInput("monthrange", "Month", 1, 12,value = 6, step =1,ticks = FALSE,sep=""),
      
      sliderInput("yearrange", "Year", 1989, 2018, value = 2000, step = 1,ticks = FALSE,sep="")    ),
    
    mainPanel(leafletOutput("mymap"))
  )
  
  
  #leafletOutput("mymap"),
  #p(),
  #sliderInput("monthrange", "Month", 1, 12,value = 6, step =1,ticks = FALSE,sep=""),
  #sliderInput("yearrange", "Year", 1989, 2018, value = 2000, step = 1,ticks = FALSE,sep="")
)

###########################################################################
#### server function
###########################################################################

server <- function(input, output, session) {
  
  output$mymap <- renderLeaflet({
    
    ### pick out which month and year you want
    month = input$monthrange[1]
    year = input$yearrange[1]
    
    iMonth=which(monthList == month & yearList == year)
    
    #### color counties by closest cell value of z
    contig_counties$SPI= zplot_ordered[,iMonth]

    #### create leaflet
    bins <- c(-3,-2,-1.5,-1,1,1.5,2,3)
    pal <- colorBin("RdYlBu", domain = contig_counties$SPI, bins = bins)
    
    leaflet(contig_counties) %>%
      addTiles() %>%
      setView(-95, 38, zoom = 4) %>%
      addPolygons(label = ~paste(NAME,": ",round(SPI,2),sep=""),
                  fillColor = ~pal(SPI),
                  color = "#444444", weight = 1,
                  smoothFactor = 0,
                  fillOpacity = 0.5,
                  highlightOptions = highlightOptions(color = "white", weight = 2,
                                                      bringToFront = TRUE)) %>%
      addLegend(pal = pal,
                values = ~SPI,
                opacity = 0.7,
                title =  paste(as.character(month),"/",as.character(year),' SPI',sep=""),
                position = "bottomright") 
    
    #### if above map has stroke=FALSE  => no lines between counties
    #### if no stroke=FALSE and color = "#444444", weight = 1, => lines between counties
    
  })
}

shinyApp(ui, server)
