\name{lpeAdj}
\alias{lpeAdj}

\title{
  High level lpeAdj function that executes the adjusted local pooled
  error significance test. If more control over parameters is needed
  then see documentation for calculateLpeAdj.
}

\description{
  Applies the LPE algorithm with two additional options.
  The first is that the original LPE method sets all variances below the
  max variance in
   the ordered distribution of variances to the maximum variance. in
   LPEadj this option is turned off by default.  The second option is to
   use a variance adjustment based on sample size rather than pi/2.
   By default the LPEadj uses the sample size based variance
   adjustment.  It is recommended to keep both of these options to the default.
}

\usage{
  lpeAdj(dat, labels=NULL, doMax=FALSE, doAdj=TRUE, q=.01)
}

\arguments{
 \item{dat}{Replicated data of experiment containing two groups (as matrix 
			or data-frame)}.
  	  
 \item{labels}{vector of group labels that correspond to the columns of
   dat. eg. labels=c(0,0,0,1,1,1) describes two groups with three
   replicates each}. 
 \item{doMax}{boolean: if T then all variances below the max variance in
   the ordered distribution of variances are set to the maximum
   variance. It is recommended to use the default value of False.}.
 \item{doAdj}{If T then run LPE with using variance adjustment value
   based on number of replicates (hardcoded in adjValues) rather than pi/2.}.
 \item{q}{ is the quantile width; q=0.01 corresponds to 100 quantiles
          i.e. percentiles. Bins/quantiles have equal number of genes and
          are split according to the average intensity A.}
 }

\details{
  The LPE test statistic numerator is the difference in medians between the
  two experimental conditions. The test statistic denominator is the combined pooled standard error for the two experimental conditions obtained by looking up the var.M from each baseOlig.error variance function. The conversion to p-values is based on the Gaussian distribution for difference if order statistics (medians). 
}

\value{
  Data frame including x, median of x, y, median of y, median difference
  of (x,y), pooled standard deviation of difference, LPE p-value,
  outlier flag, probability of an outlier within x or y, .
}
\author{Carl Murie \email{carl.murie@mcgill.ca}, 
        Nitin Jain \email{nitin.jain@pfizer.com} }


\references{
	J.K. Lee and M.O.Connell(2003). \emph{An S-Plus library for the analysis of differential expression}. In The Analysis of Gene Expression Data: Methods and Software. Edited by G. Parmigiani, ES Garrett, RA Irizarry ad SL Zegar. Springer, NewYork.

	Jain et. al. (2003) \emph{Local pooled error test for identifying
      differentially expressed genes with a small number of replicated microarrays}, Bioinformatics, 1945-1951.

Jain et. al. (2005) \emph{Rank-invariant resampling based estimation of false discovery rate for analysis of small sample microarray data}, BMC Bioinformatics, Vol 6, 187.

}

\examples{
   
 # Creating a null dataset (two groups with three
 # replicates each)
 dat <- matrix(rnorm(6000), ncol=6)
 
 # Applying LPE
 lpe.result <- lpeAdj(dat, labels=c(0,0,0,1,1,1))
  
}

\keyword{methods} % from KEYWORDS.db
