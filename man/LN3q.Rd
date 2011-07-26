% File src/library/LN3GV/man/LN3q.Rd
     \name{LN3q}
     \alias{LN3q}
     \title{Compute the expectation phase of the EM-algorithm for the LN3 model.}
     \description{Analyze microarray data using the LN3 methods detailed in Lund and Nettleton, 2011.}
     \usage{
    LN3q(parm,dat,pmat,prob,patmat)
     }
     \arguments{
	\item{parm}{Vector provides values for model parameters vg, vt, and mu.}
      \item{dat}{microarray data matrix of normalized intensities from microarray experiment.  Each row contains observations from a single gene. Each column contains observations from a single experimental unit.}
	\item{pmat}{ Matrix of posterior probabilities for expression pattern given data for each gene}
	\item{prob}{Provides prior probabilities for expression patterns}
      \item{patmat}{List containing a matrix of 1's and 0's for each expression pattern indicating which columns are equivalently expressed under the expression pattern.}
     }

\value{
The returned arguments are used within the \code{analyze} function to obtain posterior probabilities for expression patterns for each gene and model parameter estimates.
}


\author{Steve Lund \email{lunds@iastate.edu}}






