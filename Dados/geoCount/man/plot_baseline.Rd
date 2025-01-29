\name{plot_baseline}
\alias{plot_baseline}

\title{Plot Baseline Samples}

\usage{
plot_baseline(d.samples, dist.name)
}

\arguments{
  \item{d.samples}{the baseline samples}
  \item{dist.name}{the name of distance}
}

\description{
  This function plots the baseline samples.
  }



\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
plot_baseline(d.base[,1], colnames(d.base)[1])
plot_baseline(d.base[,2], colnames(d.base)[2])
plot_baseline(d.base[,3], colnames(d.base)[3])      }
}

\seealso{
\code{\link{e2dist}}, \code{\link{baseline.dist}} \code{\link{d.base}}.
}

\keyword{Transformed Residual Checking}
