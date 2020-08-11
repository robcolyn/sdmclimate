#Data inspection####

#When loading data from a database, use: 
#library(RODBC) or other appropriate database access packages
ibis <- read.csv("E://GisModelling//SDM/Ibis_SouthernBald.csv", header=TRUE) #data source/file
str(ibis)

#remove rows with missing data
library(tidyr)
#to check how many data gaps are present:
ibis.row.na <- apply(ibis, 1, function(x){any(is.na(x))})
sum(ibis.row.na)
#Drop NAs
ibis_fil1 <- ibis %>% drop_na()

#Plot presence data to inspect the extent od presence data.
#Spatial representation is important, i.e. does the spatial
#spread cover the expected range of the species? Additionally, 
#is there significant clustering? If so, spatially rarefy the 
#data.
library(raster)
#Get country codes for background country border and plot points
getData('ISO3')
countries <- getData('GADM', country='ZAF', level=1) + getData('GADM', country='SWZ', level=1) + getData('GADM', country='LSO', level=1)
map <- plot(countries, xlim=c(16,33), ylim=c(-35,-22), axes=TRUE, col="light green")
map + box() + points(ibis_fil1$Longitude, ibis_fil1$Latitude, col='black', pch=20, cex=0.75)
##alternate approach to quickly mapping data (although country borders are poor resolution):
library(maptools)
data(wrld_simpl)
map <- plot(data, xlim=c(16,33), ylim=c(-35,-22), axes=TRUE, col="light yellow")
map + box() + points(ibis_fil1$Longitude, ibis_fil1$Latitude, col='black', pch=20, cex=0.75)

#Corss-check no localities are in incorrect counteries or in the ocean
library(sp)
library(rgeos)
ibis_fil1 <- ibis_fil1[,2:3]
coordinates(ibis_fil1) <- ~Longitude+Latitude
crs(ibis_fil1) <- crs(wrld_simpl)
class(ibis_fil1)
class(wrld_simpl)
check <- over(ibis_fil1, wrld_simpl)
country <- check$NAME
outliers <- which(is.na(check))
outliers
#Once data is cleaned sufficently, convert the spatial points data 
#frame back to DF
presence_ibis <- as.data.frame(ibis_fil1)

##Generate pseudo-absence or background points####
library(dismo)
#Source directory with predictor variables, i.e. climate variables
rlist <- list.files("E://GisModelling/ClimateChange_input/Current_1970-2000/tiffs/", pattern="tif$", full.names=TRUE)
#If only bioclimatic variables are being used, this can be sourced directly from
#the Dismo package

#use the first file to create a mask, set seed for replication, generate random
#points across the enitre raster extent and/or alternatively within moving windows of 
#a user specified distance surrounding presence points
mask <- raster(rlist[1])
set.seed(1989)
backgroundpnts <- randomPoints(mask, 500) # set number of background points
backgroundpnts3 <- randomPoints(mask, 10000) #ensemble presence-only default is 10000
par(mfrow=c(1,2))
plot(!is.na(mask), legend=FALSE)
points(backgroundpnts, cex=0.5)
#Alternatively, you can constrain the background points to a specific area 
#of your raster extent (i.e. extent of occurence[EOO])
area <- extent(16, 33, -35, -22) #refine EOO lat/long here 
backgroundpnts2 <- randomPoints(mask, 50, ext=area)
plot(!is.na(mask), legend=FALSE)
plot(area, add=TRUE, col='red')
points(backgroundpnts2, cex=0.5)
#Lastly, you can further refine the locality of background selection 
#to areas within a moving window across the extent of your presence data:
projection(ibis_fil1) <- CRS('+proj=longlat +datum=WGS84')
#create circles with a radius that is meaningful to the species in 
#question, 150 km used here as it represents the dispersal ability. 
windows <- circles(ibis_fil1, d=150000, lonlat=TRUE)
circles <- polygons(windows)
#Random smapling within circles. For stratified or other supported sampling types
#change "type". 
sample1 <- spsample(circles, 250, type='random', iter=25)
cells <- cellFromXY(mask, sample1)
length(cells)
cells <- unique(cells)
length(cells)
xy <- xyFromCell(mask, cells)
plot(circles, axes=TRUE)
points(xy, cex=0.75, pch=20, col='blue') + points(ibis_fil1$Longitude, ibis_fil1$Latitude, col='black', pch=20, cex=0.75)

#Convert from spatial points DF back to DF
background_ibis <- as.data.frame(backgroundpnts3) #change "backgroundpnts3"




