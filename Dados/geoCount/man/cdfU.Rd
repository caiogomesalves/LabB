\name{cdfU}
\alias{cdfU}

\title{Approximate the CDF Value from Reference Samples}

\usage{
cdfU(Y.obs, Y.rep, discrete = FALSE)
}

\arguments{
  \item{Y.obs}{a vector which indicates the observed data set}
  \item{Y.rep}{a matrix which indicates the reference data sets}
  \item{discrete}{a logical value which indicates if the variable is discrete}
}

\description{
  This function approximates the CDF value for the observed data by using reference data.}


\value{
 A vector of CDF values.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Y.obs <- 11:20
res <- matrix(0, 10, 50)
for(i in 1:50){
Y.rep <- matrix(rpois(10*5000, 15), 10, )
res[, i] <- cdfU(Y.obs, Y.rep)
}
matplot(t(res), type="l")
abline(h = ppois(11:20, 15))
}
}

\seealso{
\code{\link{tranR}}
}

\keyword{Transformed Residual Checking}
