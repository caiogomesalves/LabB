## R code from Section 1.3, "Maps and geodesics in R," pp.16--19

library(maps)

mn.map = map(database="county", region="minnesota")

## If we do not want county boundaries, then

mn.map = map(database="state", region="minnesota")

## Importing shapefiles
library(maptools)

minnesota.shp <- readShapePoly("minnesota.shp",
       proj4string=CRS("+proj=longlat"))

plot(minnesota.shp)

## Map of Canada
library(mapdata)

map("worldHires", "Canada", xlim=c(-141,-53), ylim=c(40,85), col="gray90", fill=TRUE)


## Mapping Influenza A (H1N1) cases in Minnesota
newdata <- read.csv("http://www.biostat.umn.edu/~sudiptob/JSM_2014/R_Programs/newdata.csv")

minnesota.shp@data <- merge(minnesota.shp@data, newdata, by="NAME")

library(RColorBrewer)
library(classInt)

CASES <- minnesota.shp@data$cases
POP1999 <- minnesota.shp@data$POP1999
var <- (CASES/POP1999)*100
nclr <- 7
plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(var, nclr, style="equal", dataPrecision=4)
colcode <- findColours(class, plotclr)
plot(minnesota.shp)
plot(minnesota.shp, col=colcode, add=T)
title(main="Influenza A (H1N1) - Rate per 100")
legend("bottomright", legend=names(attr(colcode, "table")),
    fill=attr(colcode, "palette"), cex=0.6, bty="n")


## Map projections in R
library(mapproj)
sinusoidal.proj = map(database= "world", ylim=c(45,90), xlim=c(-160,-50), col="grey80", fill=TRUE, projection="sinusoidal", plot=FALSE)
map(sinusoidal.proj)

## For distance computations using map projections
coloradoST <- read.table("http://www.biostat.umn.edu/~sudiptob/JSM_2014/R_Programs/ColoradoS-T.dat", header=TRUE)
LON <- coloradoST$Longitude
LAT <- coloradoST$Latitude

xy.sinusoidal <- mapproject(LON, LAT, projection="sinusoidal")
xy.coords = cbind(xy.sinusoidal$x, xy.sinusoidal$y)

## Comparing Geodesic and Euclidean distances in the Colorado dataset
library(fields)
##scale planar coordinates by radius of earth
radius.of.earth = 6371 ##kilometers or 3963.17miles
XY.sinusoidal = radius.of.earth*(cbind(xy.sinusoidal$x, xy.sinusoidal$y))
euclidean.dist = rdist(XY.sinusoidal)
XY.sphere = cbind(LON, LAT) ##First column must be longitude, second is latitude
spherical.dist = rdist.earth(XY.sphere, miles=FALSE)


## Mapping with RgoogleMaps
library(RgoogleMaps)

long <- LON
lat <- LAT
MyMap <- GetMap.bbox(lonR = range(long), latR = range(lat),
                     size=c(640,640), maptype = "hybrid")

PlotOnStaticMap(MyMap)

#convert to long/lat to the same coordinate as in MyMap
convert_points <- LatLon2XY.centered(MyMap, lat, long)
points(convert_points$newX, convert_points$newY, col = 'red', pch=19)


##Converting Long/Lat to utm#
library(sp)
library(rgdal)

SP_longlat <- SpatialPoints(coords = cbind(coloradoST$Longitude,coloradoST$Latitude), proj4string = CRS("+proj=longlat +ellps=WGS84"))
SP_utm <- spTransform(SP_longlat, CRS("+proj=utm +zone=18 +datum=WGS84"))
plot(SP_utm)



