\name{cutChain}
\alias{cutChain}

\title{Modify Markov Chains with Burn-in and Thining}

\usage{
cutChain(res, chain.ind=2:4, burnin, thinning)
}

\arguments{
  \item{res}{a list with elements containing the posterior samples of latent variables and parameters; usually the output from \code{\link{runMCMC}} }
  \item{chain.ind}{the index of elements in "res" that will be modified}
  \item{burnin}{the number for burn-in}
  \item{thinning}{the number for thining}
}

\description{
  This function takes the results from \code{\link{runMCMC}} and modifies the chains of posterior samples for burn-in and thinning.}


\value{
 A list with elements containing the modified posterior samples.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
res <- runMCMC(Y, L=0, loc=loc, MCMCinput = input )
res.m <- cutChain(res, chain.ind=1:4, burnin=100, thining=10)
}
}

\seealso{
   \code{\link{runMCMC}}, \code{\link{mixChain}}.
}

\keyword{MCMC}
