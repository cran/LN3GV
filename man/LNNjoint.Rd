% File src/library/LN3GV/man/LNNjoint.Rd
     \name{LNNjoint}
     \alias{LNNjoint}
     \title{Fit a hierarchical model to matrix of normalized microarray data}
     \description{Analyze microarray data using the LNN methods detailed in Lund and Nettleton, 2011.}
     \usage{
    LNNjoint(dat,prob,parm,patmat)
     }
     \arguments{
      \item{dat}{microarray data matrix of normalized intensities from microarray experiment.  Each row contains observations from a single gene. 
Each column contains observations from a single experimental unit.}
	\item{prob}{Provides prior probabilities for expression patterns}
	\item{parm}{Vector provides values for model parameters vt, vr, and mu.}
      \item{patmat}{List containing a matrix of 1's and 0's for each expression pattern indicating which columns are equivalently expressed under the expression pattern.}

     }

\value{
The returned arguments are used within the \code{analyze} function to obtain posterior probabilities for expression patterns for each gene and model parameter estimates.

list containing:
	\item{"Joint"}{Matrix providing joint density of expression pattern and data for each gene }
	\item{"Marginal"}{vector marginal density for data for each gene.  (This is the row sum of `Joint'.)}
	\item{"Posterior"}{Posterior probabilities for expression pattern given data (by gene)}
	}

\author{Steve Lund \email{lunds@iastate.edu}}






