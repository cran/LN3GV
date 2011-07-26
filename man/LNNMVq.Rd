% File src/library/LN3GV/man/LNNMVq.Rd
     \name{LNNMVq}
     \alias{LNNMVq}
     \title{Compute the expectation phase of the EM-algorithm for the LNNMV model.}
     \description{Analyze microarray data using the LNNMV* methods detailed in Lund and Nettleton, 2011.}
     \usage{
    LNNMVq(parm,dat,pmat,prob,vr,pats)
     }
     \arguments{
	\item{parm}{Vector provides values for model parameters vg, vt, and mu.}
      \item{dat}{microarray data matrix of normalized intensities from microarray experiment.  Each row contains observations from a single gene. Each column contains observations from a single experimental unit.}
	\item{pmat}{ Matrix of posterior probabilities for expression pattern given data for each gene}
	\item{prob}{Provides prior probabilities for expression patterns}
	\item{vr}{Vector providing gene-specific error variances}
      \item{pats}{matrix in which each row describes a unique expression pattern. Must have same number of columns as \code{dat}.}
     }

\value{
The returned arguments are used within the \code{analyze} function to obtain posterior probabilities for expression patterns for each gene and model parameter estimates.
}

\author{Steve Lund \email{lunds@iastate.edu}}






