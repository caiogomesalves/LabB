\name{d.base}
\alias{d.base}

\title{Data Set of Baseline Samples}

\usage{
data(Dbase_n100N5000)
}


\description{
  This data set contains baseline samples for 100 residuals with 5000 iterations.}

\details{
  A data.frame "d.base" with 5000 observations and 3 variables will be loaded.
} 


\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{

data(Dbase_n100N5000)
str(d.base) 
plot_baseline(d.base[,1], colnames(d.base)[1])
}
}

\seealso{
\code{\link{baseline.dist}}, \code{\link{plot_baseline}}, \code{\link{pOne}}.
}

\keyword{Transformed Residual Checking}
