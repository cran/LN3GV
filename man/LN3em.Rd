% File src/library/LN3GV/man/LN3em.Rd
     \name{LN3em}
     \alias{LN3em}
     \title{Fit the LN3 model to matrix of normalized microarray data}
     \description{Analyze microarray data using the LN3 methods detailed in Lund and Nettleton, 2011.}
     \usage{
    LN3em(dat,pats,MaxIt,parm.conv,prob.conv,
trt,parm.init=NULL,prob.init=NULL)
     }
     \arguments{
      \item{dat}{microarray data matrix of normalized intensities from microarray experiment.  Each row contains observations from a single gene. 
Each column contains observations from a single experimental unit.}
      \item{pats}{matrix in which each row describes a unique expression pattern. Must have same number of columns as \code{dat}.  For each pattern, experimental units sharing a common integer are assumed to be equivalently expressed.}
	\item{MaxIt}{maximum number of EM algorithm iterations used to estimate prior probabilities for expression patterns and model parameters.}
	\item{parm.conv}{Until iteration \code{MaxIt}, EM algorithm will not stop while the difference between consecutive iterations in estimates for any parameter is greater than \code{parm.cov} times the updated parameter estimate. May be specified as single number or a vector with length equal to number of model parameters (see details).}
	\item{prob.conv}{Until iteration \code{MaxIt}, EM algorithm will not stop while the difference between consecutive iterations in prior probability estimates for any pattern is greater than \code{prob.conv}. May be specified as single number or a vector with length equal to \code{nrow(pats)}.}
	\item{trt}{Vector of integers assigning columns to groups.  Columns sharing a common assigned integer are assumed to be replicates of one another and therefore equivalently expressed with each other for all genes.}
	\item{parm.init}{Optional. Vector providing initial estimates for model parameters vg, vt, vr, and mu.}
	\item{prob.init}{Optional. Provides initial estimates of prior probabilities for expression patterns}

     }

\value{
The returned arguments are used within the \code{analyze} function to obtain posterior probabilities for expression patterns for each gene and model parameter estimates.

list containing:
	\item{"Pattern Probs"}{Vector of prior probabilities for expression patterns}
	\item{"Parameter Estimates"}{vector providing model parameters vg, vt, vr, and mu}

}

\author{Steve Lund \email{lunds@iastate.edu}}






