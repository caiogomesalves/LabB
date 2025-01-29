\name{plotData}
\alias{plotData}

\title{Plot Geostatistical Data}

\usage{
plotData(Y=NULL, loc=NULL, Yp=NULL, locp=NULL, Yt=NULL, loct=NULL, bdry=NULL, 
         cols=1:3, pchs=1:3, size=c(0.3, 2.7), ...)
}

\arguments{
  \item{Y, Yp, Yt}{the vector of response variables}
  \item{loc, locp, loct}{\eqn{n \times 2} matrix that indicates the coordinates of locations}
  \item{bdry}{a list containing the coordinates of boundaries}
  \item{cols}{the colors used for different sets of response variables}
  \item{pchs}{the shapes used for different sets of response variables}
  \item{size}{the minimum and maximum of the sizes}
  \item{...}{other parameters that control the plotting}
}

\description{
  This function plots geostatistical data for up to three data sets.
}


\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
loc1 <- locGrid(1, 1, 10, 10)
loc2 <- locCircle(0.72, 60)
loc3 <- locSquad(0.38, 10)
loc <- rbind(as.matrix(loc1) , loc2, loc3); plot(loc)
dat <- simData(loc, cov.par = c(1, 0.2, 1))
Y <- dat$data
plotData(Y[1:nrow(loc1)], loc1, 
         Y[(nrow(loc1)+1):(nrow(loc1)+nrow(loc2))], loc2, 
         Y[(length(Y)-nrow(loc3)+1):length(Y)], loc3, 
         xlab="x", ylab="y", pchs = c(1, 16, 15)
         )
# plot boundaries
data(TexasCounty_boundary)
plotData(bdry = TexasCounty.boundary)
# plot data with the boundary
data(Rongelap)
str(Rongelap)
plotData(bdry = Rongelap$borders, Y = Rongelap$data, loc = Rongelap$coords)
}
}

\seealso{
   \code{\link{locCircle}},
   \code{\link{locGrid}},
   \code{\link{locSquad}}, \code{\link{simData}}.
}

\keyword{Data}
