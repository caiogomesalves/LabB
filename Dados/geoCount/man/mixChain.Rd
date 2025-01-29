\name{mixChain}
\alias{mixChain}

\title{Mix Parallel Markov Chains}

\usage{
mixChain(res.m.prl)
}

\arguments{
  \item{res.m.prl}{a list with each element containing the result of posterior samples from one CPU; the elements shoul only contain the Markov chains of posterior samples (while "AccRate" is eliminated when using \code{\link{cutChain}}) }
}

\description{
  This function mix parallel chains into one chain.}


\value{
 A list with elements containing the mixed posterior samples.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
res.prl <- runMCMC.multiChain(Y, L=0, loc=loc, X=loc, 
            MCMCinput = input, n.chn = 5)
res.m.prl <- lapply(res.prl, cutChain, chain.ind=1:4, burnin=200, thining=20)
res.mix <- mixChain(res.m.prl)
}
}

\seealso{
   \code{\link{cutChain}}.
}

\keyword{MCMC}
