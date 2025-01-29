\name{plot_etran}
\alias{plot_etran}

\title{Plot Transformed Residuals}

\usage{
plot_etran(e.tran, fig = 1:4)
}

\arguments{
  \item{e.tran}{a vector which indicates the transformed residuals}
  \item{fig}{a vector which indicates which types to plot: 1 indicates scatter plot, 2 indicates QQ-plot, 3 indicates density plot, and 4 indicates relative density plot (with standard normal distribution served as the base)}
}

\description{
  This function plots transformed residuals in different types.}

\details{
  \code{density} function is used to compute the empirical density.
  
  \code{reldist} function in \code{{reldist}} is used to compute the relative density.
} 


\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
require(reldist)
plot_etran(rnorm(200), fig=c(1,4))   }
}

\seealso{
\code{\link{tranR}}
}

\keyword{Transformed Residual Checking}
