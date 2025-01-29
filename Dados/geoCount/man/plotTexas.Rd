\name{plotTexas}
\alias{plotTexas}

\title{Plot Texas Counties}

\usage{
plotTexas(TexasCounty.boundary, ind.col = sample(2:5, 254, replace=T))
}

\arguments{
  \item{TexasCounty.boundary}{a list containing the boundary data of Texas counties, see example for details}
  \item{ind.col}{a vector used for plotting parameter "col" to indicate the colors of counties}  
}

\description{
  This function plot all the 254 Texas counties.}

\details{
  This function uses \code{polygon} function to draw the boundaries.
} 


\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
data(TexasCounty_boundary)
plotTexas(TexasCounty.boundary)
}
}


\keyword{Data}
