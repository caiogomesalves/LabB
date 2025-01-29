\name{Weed}
\docType{data}
\alias{Weed}

\title{Weed Data}

\usage{Weed}


\description{
  This data set contains the Weed data set.}

\details{
  The Weed data were collected at the Bjertorp farm in the south-west of Sweden. Weed counts
of non-crop plants were observed at different locations, and camera recorded images were
used to estimate the counts with the help of certain image analysis algorithm. Guillot, G., Loren, N., and Rudemo, M. (2009). Spatial prediction of weed intensities from exact count data and image-based estimates. Journal of Applied Statistics, 58, Part 4, 525-542.
} 

\value{
  A data.frame with 100 observations and 4 variables.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
data(Weed)
str(Weed)
plotData(Weed[,3], Weed[,1:2], Weed[,4], Weed[,1:2], 
         xlab="Eastings", ylab="Northings")
}
}

\seealso{
\code{\link{plotData}}.
}

\keyword{Data}
