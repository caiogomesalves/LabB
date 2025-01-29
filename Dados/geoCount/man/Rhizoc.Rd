\name{Rhizoc}
\docType{data}
\alias{Rhizoc}

\title{Rhizoc Data}

\usage{Rhizoc}


\description{
  This data set contains the Rhizoc data set.}

\details{
This dataset appeared Zhang (2002) and was modeled 
using the Binomial Logit-normal spatial model.
It consists of counts of the root disease Rhizoctonia root rot present 
in barey plants collected at 100 locations in the Cunningham farm in the
north-west of the U.S. 
The sampling consisted on pulling from the ground 15 plants at each location.
This dataset is a 100 by 4 matrix with rows \eqn{(x_i , y_i , t_i)}, 
where \eqn{x_i} are the coordinates of the i\eqn{^{th}} sampling location, 
\eqn{y_i} is the total number of infected crown roots in the pulled plants at 
\eqn{x_i}, and \eqn{t_i} is the total number of crown roots in the pulled plants at \eqn{x_i}; see  Zhang (2002) for further details.
} 

\value{
  A data.frame with 100 observations and 5 variables.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
data(Rhizoc)
str(Rhizoc)
plotData(Rhizoc[,3], Rhizoc[,1:2], Rhizoc[,4], Rhizoc[,1:2],
         xlab="Eastings", ylab="Northings")
}
}

\seealso{
\code{\link{plotData}}.
}

\keyword{Data}
