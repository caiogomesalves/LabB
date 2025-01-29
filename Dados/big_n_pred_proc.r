
## R programs illustrating multivariate predictive process in BCG Second Edition, pages 405--410

library(spBayes)
library(fields)
library(cluster)
library(rgdal)
library(maptools)

## We begin by removing non-forest inventory plots, converting biomass measurements from kilograms per hectare to the log of metric tons per hectare, and taking a look at plot locations across the forest

data(BEF.dat)

BEF.dat <- BEF.dat[BEF.dat$ALLBIO02_KGH>0,]
bio <- BEF.dat$ALLBIO02_KGH*0.001

bio <- BEF.dat$ALLBIO02_KGH*0.001;
log.bio <- as.matrix(log(0.001*BEF.dat[, c("BOLE02_KGH","BRANCH02_KGH","FOLIAGE02_KGH")]))
colnames(log.bio) <- c("log.bole.mt", "log.branch.mt","log.foliage.mt")
coords <- as.matrix(BEF.dat[,c("XUTM","YUTM")])
plot(coords, pch=19, cex=0.5, xlab="Easting (m)", ylab="Northing (m)")

## To get an idea of the covariances among the three outcomes
cov(log.bio)

## Examples of setting up knots using k-means, k-medoids and cover.design in R
m <- 50
km.knots <- kmeans(coords, m)$centers
cl.knots <- clara(coords, m)$medoids
cd.knots <- cover.design(coords, nd=m)$design

## Plot the knots to see how well they cover the domain
plot(coords, pch=19, cex=0.5, xlab="Easting (m)", ylab="Northing (m)")
points(km.knots, pch=5, cex=1, col="blue")
points(cl.knots, pch=6, cex=1, col="green")
points(cd.knots, pch=7, cex=1, col="red")
legend("bottomright", cex=1, pch=c(19,5,6,7), bty="n", col=c("black","blue","green","red"), legend=c("observations","kmeans","clara", "cover.design"))

## Create a model "list" for the multivariate spatial regression and run spMvLM using the modified predictive process with knots from k-means
q <- 3
n.samples <- 2000
n.ltr <- q*(q+1)/2
model <- list(log.bio[,"log.bole.mt"] ~ ELEV + SLOPE + SUM_02_TC1 + SUM_02_TC2 + SUM_02_TC3, log.bio[,"log.branch.mt"]~ ELEV + SLOPE + SUM_02_TC1 + SUM_02_TC2 + SUM_02_TC3, log.bio[,"log.foliage.mt"] ~ ELEV + SLOPE + SUM_02_TC1 + SUM_02_TC2 + SUM_02_TC3)

bef.spMvLM <- spMvLM(model, coords=coords, knots=km.knots, data=BEF.dat, starting=list("phi"=rep(3/500,q), "A"=rep(0.05,n.ltr), "Psi"=rep(0.05,q)), tuning=list("phi"=rep(0.3,q), "A"=rep(0.0001,n.ltr), "Psi"=rep(0.001,q)), priors=list("phi.Unif"=list(rep(3/2500,q), rep(3/100,q)), modified.pp=TRUE, "K.IW"=list(q+1, diag(0.1,q)), "Psi.IG"=list(rep(2,q), rep(0.05,q))), cov.model="exponential", n.samples=n.samples, verbose=TRUE, n.report=100)

## Set burnin. Recover regression parameters and spatial random effects
burn.in <- floor(0.75*n.samples)
bef.spMvLM <- spRecover(bef.spMvLM, start=burn.in)
round(summary(mcmc(cbind(bef.spMvLM$p.beta.recover.samples, bef.spMvLM$p.theta.recover.samples)))$quantiles[,c(1,3,5)],3)

## Obtain fitted values and plot them
fitted <- spPredict(bef.spMvLM, start=burn.in, thin=10, pred.covars=bef.spMvLM$X, pred.coords=bef.spMvLM$coords)
y.hat <- rowMeans(fitted$p.y.predictive.samples)
bole <- y.hat[seq(1,length(y.hat),q)]
branch <- y.hat[seq(2,length(y.hat),q)]
foliage <- y.hat[seq(3,length(y.hat),q)]

res <- 100
par(mfrow=c(2,2))
surf <- mba.surf(cbind(coords,bole), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Bole fitted values")
points(coords)
surf <- mba.surf(cbind(coords,branch), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Branch fitted values")
points(coords)
surf <- mba.surf(cbind(coords,foliage), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Foliage fitted values")
points(coords)


## Spatial predictions

x.range <- range(coords[,1])
y.range <- range(coords[,2])
BEF.shp <- readShapePoly("BEF-data/BEF_bound.shp")
BEF.poly <- as.matrix(BEF.shp@polygons[[1]]@Polygons[[1]]@coords)
BEF.grids <- readGDAL("BEF-data/dem_slope_lolosptc_clip_60.img")

pred.covars <- cbind(BEF.grids[["band1"]],BEF.grids[["band2"]], BEF.grids[["band3"]],BEF.grids[["band4"]], BEF.grids[["band5"]])
pred.covars <- cbind(rep(1, nrow(pred.covars)), pred.covars)
pred.coords <- SpatialPoints(BEF.grids)@coords
pred.covars <- pred.covars[pointsInPoly(BEF.poly, pred.coords),]
pred.coords <- pred.coords[pointsInPoly(BEF.poly, pred.coords),]
pred.X <- mkMvX(list(pred.covars, pred.covars, pred.covars))
bef.bio.pred <- spPredict(bef.spMvLM, start=burn.in, thin=5, pred.coords=pred.coords, pred.covars=pred.X)

y.pred.mu <- apply(bef.bio.pred$p.y.predictive.samples, 1, mean)
y.pred.sd <- apply(bef.bio.pred$p.y.predictive.samples, 1, sd)
bole.pred.mu <- y.pred.mu[seq(1,length(y.pred.mu),q)]
branch.pred.mu <- y.pred.mu[seq(2,length(y.pred.mu),q)]
foliage.pred.mu <- y.pred.mu[seq(3,length(y.pred.mu),q)]
bole.pred.sd <- y.pred.sd[seq(1,length(y.pred.sd),q)]
branch.pred.sd <- y.pred.sd[seq(2,length(y.pred.sd),q)]
foliage.pred.sd <- y.pred.sd[seq(3,length(y.pred.sd),q)]

bio.pred.grid <- as.data.frame(list(x=pred.coords[,1], y=pred.coords[,2], bole.mu=bole.pred.mu, branch.mu=branch.pred.mu, foliage.mu=foliage.pred.mu, bole.sd=bole.pred.sd, branch.sd=branch.pred.sd, foliage.sd=foliage.pred.sd))
coordinates(bio.pred.grid) <- c("x", "y") ## to SpatialPointsDataFrame
gridded(bio.pred.grid) <- TRUE ## to SpatialGridDataFrame

toImage <- function(x){as.image.SpatialGridDataFrame(x)}
res <- 100
par(mfrow=c(3,2))
surf <- mba.surf(cbind(coords,log.bio[,"log.bole.mt"]), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(surf[["z"]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", main="Observed bole biomass")
image.plot(toImage(bio.pred.grid["bole.mu"]), zlim=z.lim, xaxs = "r", yaxs = "r", main="Predicted bole biomass")
surf <- mba.surf(cbind(coords,log.bio[,"log.branch.mt"]), no.X=res, no.Y=res,extend=FALSE)$xyz.est
z.lim <- range(surf[["z"]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", main="Observed branch biomass")
image.plot(toImage(bio.pred.grid["branch.mu"]), zlim=z.lim, xaxs = "r", yaxs = "r", main="Predicted branch biomass")
surf <- mba.surf(cbind(coords,log.bio[,"log.foliage.mt"]), no.X=res, no.Y=res, extend=FALSE)$xyz.est
z.lim <- range(surf[["z"]], na.rm=TRUE)
image.plot(surf, xaxs = "r", yaxs = "r", main="Observed foliage biomass")
image.plot(toImage(bio.pred.grid["foliage.mu"]), zlim=z.lim, xaxs = "r", yaxs = "r", main="Predicted foliage biomass")





