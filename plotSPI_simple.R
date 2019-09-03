####  simpleSPIplot.R
####
####  very simple script to plot the results of the SPI calculation
####

library(mapdata)
library(fields)

load("SPI_out_3month.RData")

### pick out which estimate you want to use z, z_emp, z_fit
zplot <- spi_results$z_emp

### pick out which month and year you want
month = 6
year = 2006

iMonth=which(spi_results$monthList == month & spi_results$yearList == year)
image.plot(Lon,Lat,zplot[,,iMonth],xlim=range(-125,-66.5),ylim=range(24.5,50),zlim=range(-3,3),col=rainbow(7),breaks=c(-3,-2,-1.5,-1,1,2,2.5,3))
map("county",add=TRUE)
