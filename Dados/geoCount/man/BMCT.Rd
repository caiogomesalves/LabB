\name{BMCT}
\alias{BMCT}

\title{Perform Bayesian Model Checking}

\usage{
BMCT(Y.obs, Y.rep, funcT, ifplot = FALSE)
}

\arguments{
  \item{Y.obs}{a vector which indicates the observed data set}
  \item{Y.rep}{a matrix which indicates the reference data sets}
  \item{funcT}{a function which defines the dignostic statistic}
  \item{ifplot}{a logical value which indicates whether plot the dignostic statistics}
}

\description{
  This function conducts Bayesian model checking by comparing observed and reference data sets and reveals the result via "p-value" and "RPS" (as well as the plot).
  }


\value{
 A vector of p-value and RPS.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Yrep.eb <- repYeb(N.sim=2000, loc, L, res.m, est = "mode")
funcT <- function(Y){ max(Y)-min(Y) }
BMCT(Y, Yrep.eb, funcT, ifplot=TRUE)  }
}

\seealso{
\code{\link{repYeb}}, \code{\link{repYpost}}, \code{\link{pRPS}}, \code{\link{plot_pRPS}}
}

\keyword{Bayesian Model Checking}
