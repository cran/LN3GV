% File src/library/LN3GV/man/analyze.Rd
     \name{analyze}
     \alias{analyze}
     \title{Fit a hierarchical model to matrix of normalized microarray data}
     \description{
       Analyze microarray data using the methods detailed in Lund and Nettleton, 2011.  This is the main function of the LN3GV package.}
     \usage{
    analyze(dat,pats,method,log.scale=TRUE,MaxIt=30,
parm.conv=.005,prob.conv=.0005,parm.init=NULL,prob.init=NULL)
     }
     \arguments{
      \item{dat}{microarray data matrix of normalized intensities from microarray experiment.  Each row contains observations from a single gene. 
Each column contains observations from a single experimental unit.}
      \item{pats}{matrix in which each row describes a unique expression pattern. Must have same number of columns as \code{dat}.  For each pattern, experimental units sharing a common integer are assumed to be equivalently expressed.}
	\item{method}{method used to analyze data. Must be one of "LNNMV*", "LNNGV", "LN3", "LN3MV*", or "LN3GV"}
	\item{log.scale}{logical.  If false, function will perform analysis using \code{log(dat)}}
	\item{MaxIt}{maximum number of EM algorithm iterations used to estimate prior probabilities for expression patterns and model parameters.}
	\item{parm.conv}{Until iteration \code{MaxIt}, EM algorithm will not stop while the difference between consecutive iterations in estimates for any parameter is greater than \code{parm.cov} times the updated parameter estimate. May be specified as single number or a vector with length equal to number of model parameters (see details).}
	\item{prob.conv}{Until iteration \code{MaxIt}, EM algorithm will not stop while the difference between consecutive iterations in prior probability estimates for any pattern is greater than \code{prob.conv}. May be specified as single number or a vector with length equal to \code{nrow(pats)}.}
	\item{parm.init}{Optional. Provides initial estimates for model parameters.  See details. }
	\item{prob.init}{Optional. Provides initial estimates of prior probabilities (or mixing proportions) for expression patterns.}
     }

\details{
Order of model parameters optimized using EM algorithm for each method: 

LNNMV* and LNNGV: Treatment variance (vt), mu.  Error variance parameter values are not chosen using EM algorithm for these methods. 

LN3: Gene variance (vg), Treatment variance (vt), Error variance (vr), mu 

LN3MV* and LN3GV: Gene variance (vg), Treatment variance (vt), mu. Error variance parameter values are not chosen using EM algorithm for these methods. 


When estimating error variances, experimental units for which the corresponding columns of \code{pats} are identical are assumed to be replicates.
     }

\value{list containing:
	\item{"Pattern Probs"}{matrix the describes posterior probabilities for expression patterns for each gene}
	\item{"Parameter Estimates"}{vector providing estimated prior probabilities for expression patterns and model parameters.}
}

\author{Steve Lund \email{lunds@iastate.edu}}



     \examples{  
### Create example data set from LN3MV model.  
dat<-matrix(rnorm(500),50,10)*rgamma(50,.5,1)+rnorm(50,0,2); 

### Make first 25 example genes differentially expressed.
### Suppose 2 conditions with 4 and 6 reps, respectively.
dat[1:25,1:4]<-dat[1:25,1:4]+rnorm(25)
dat[1:25,5:10]<-dat[1:25,5:10]+rnorm(25)
dat[26:50,]<-dat[26:50,]+rnorm(25)

### Create matrix defining possible expression patterns.
### First row is pattern for equivalent expression.
### Second row is pattern for differential expression across conditions. 
pats<-rbind(rep(1,10),rep(1:2,c(4,6)))

### Analyze data using each method
analyze(dat,pats,method="LN3GV")
analyze(dat,pats,method="LN3MV*")
analyze(dat,pats,method="LN3")
analyze(dat,pats,method="LNNGV")
analyze(dat,pats,method="LNNMV*")
}
     \keyword{microarray, differential expression}



