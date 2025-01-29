\name{Rongelap}
\docType{data}
\alias{Rongelap}

\title{Data Set of Rongelap Island}

\usage{Rongelap}


\description{
  This data set contains the Rongelap data.}

\details{
  The data were collected from Rongelap Island, the principal island of Rongelap Atoll in the South Pacific, which forms part of the Marshall Islands. U.S. nuclear weapon testing program generated heavy fallout over the island in the 1950s and it has been uninhabited since 1985. Diggle, P. J., Tawn, J. A. and Moyeed, R. A. (1998). Model based geostatistics (with discussion). Applied Statistics, 47, 299-350.
} 

\value{
A list with 4 elements: 
  \item{coords}{a \eqn{157 \times 2} matrix which indicates the coordinates of 157 sampled locations}
  \item{data}{a vector of length 157 indicates the counts of photo emission for 157 sampled locations}
  \item{units.m}{a vector of length 157 indicates the time (in seconds) over which the counts was accumulated}
  \item{borders}{a matrix containing the boundary information of Rongelap island}
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
data(Rongelap)
str(Rongelap)
plotDataBD(Rongelap$borders, Rongelap$data, Rongelap$coords)
}
}

\seealso{
\code{\link{plotData}}.
}

\keyword{Data}
