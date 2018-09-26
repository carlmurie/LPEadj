\name{calculateLpeAdj}
\alias{calculateLpeAdj}

\title{
  Evaluates local pooled error significance test with user chosen
  variance adjustments.
}

\description{
  The local pooled error test attempts to reduce dependence on the
  within-gene estimates in tests for differential expression, by pooling
  error estimates within regions of similar intensity. Note that with
  the large number of genes there will be genes with low within-gene
  error estimates by chance, so that some signal-to-noise ratios will be
  large regardless of mean expression intensities and fold-change. The
  local pooled error attempts to avert this by combining within-gene
  error estimates with those of genes with similar expression
  intensity.
}

\usage{
  calculateLpeAdj(x, y, basevar.x,basevar.y, df=10, array.type="olig", 
      probe.set.name="OLIG.probe.name", trim.percent=5, adjust1=1.57, adjust2=1.57)
}

\arguments{
 \item{x}{Replicated data from first experimental condition (as matrix 
			or data-frame)}.
  	  
 \item{y}{Replicated data from second experimental condition (as matrix 
			or data-frame)}.
 \item{basevar.x}{Baseline distribution of first condition obtained from 
		 function baseOlig.error}
 \item{basevar.y}{Baseline distribution of second condition obtained from 
		 function baseOlig.error}
 \item{df}{Degrees of freedom used in fitting smooth.spline to estimates
	   of var.M for bins in A}
 \item{array.type}{Currently supports oligo arrays}
 \item{probe.set.name}{Gene IDs. By default if they are not provided
   then 1,2,3,... is assigned as GeneID}
 \item{trim.percent}{Percent of (A, var.M) estimates to trim from low 
   end of A}
 \item{adjust1}{adjustment factor of variance for first group}
 \item{adjust2}{adjustment factor of variance for second group}
 }

\details{
  The LPE test statistic numerator is the difference in medians between the
  two experimental conditions. The test statistic denominator is the combined pooled standard error for the two experimental conditions obtained by looking up the var.M from each baseOlig.error variance function. The conversion to p-values is based on the Gaussian distribution for difference if order statistics (medians). The user may select both the smoother degrees of freedom (smaller is smoother) and the trim percent to obtain a variance function to suit particular 
issues i.e. variability of genes with low expression intensity. The
default values for the adjustment of the variances of the two groups is
the asymptotically correct value of pi/2.  This value is biased at small sample
values and unbiased adjustment parameters based on sample size can be
used instead. See documentation of lpeAdj for details.
}

\value{
  Data frame including x, median of x, y, median of y, median difference of (x,y), pooled standard deviation of difference, LPE p-value, outlier flag, probability of an outlier within x or y.
}
\author{ Carl Murie \email{carl.murie@mcgill.ca},
         Nitin Jain \email{nitin.jain@pfizer.com}}


\references{
	J.K. Lee and M.O.Connell(2003). \emph{An S-Plus library for the analysis of differential expression}. In The Analysis of Gene Expression Data: Methods and Software. Edited by G. Parmigiani, ES Garrett, RA Irizarry ad SL Zegar. Springer, NewYork.

	Jain et. al. (2003) \emph{Local pooled error test for identifying
      differentially expressed genes with a small number of replicated microarrays}, Bioinformatics, 1945-1951.

Jain et. al. (2005) \emph{Rank-invariant resampling based estimation of false discovery rate for analysis of small sample microarray data}, BMC Bioinformatics, Vol 6, 187.

}

\seealso{
  \code{\link{lpeAdj}}
}

\examples{
   
 # Loading the data from the LPE library
 data(Ley)

 ADJ.VALUES <- c(1, 1, 1.34585905516761 ,1.19363228146169 ,1.436849413109
                  ,1.289652132873 ,1.47658053092781 ,1.34382984852146
                  ,1.49972130857404, 1.3835405678718)
            
 dim(Ley)
 # Gives 12488*7 
 # First column is ID.



 # Subsetting the data
 subset.Ley <- Ley[1:1000,]
  
 subset.Ley[,2:7] <- preprocess(subset.Ley[,2:7],data.type="MAS5")
  
 # Finding the baseline distribution of condition 1 and 2.
 var.1 <- adjBaseOlig.error(subset.Ley[,2:4], q=0.01, setMax1=FALSE)
 var.2 <- adjBaseOlig.error(subset.Ley[,5:7], q=0.01, setMax1=FALSE)
 
 # Applying LPE
 lpe.result <- calculateLpeAdj(subset.Ley[,2:4],subset.Ley[,5:7], var.1, var.2,
		probe.set.name=subset.Ley[,1], adjust1=ADJ.VALUES[3],
                adjust2=ADJ.VALUES[3])
  
}

\keyword{methods} % from KEYWORDS.db
