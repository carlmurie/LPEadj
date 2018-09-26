library(LPE)

#
# This file is a correction for the LPE statistical test (Jain et al, 2003).
# It does not modify the algorithm in any way apart from the two steps listed
# below.  All LPE documentation is still correct apart from what is shown
# below.
#
# The correction is in two parts.
#
# 1. In the LPE bioconductor package the LPE method sets all variances below
#    the max variance in the ordered distribution of variances to the maximum
#    variance. This step has been shown to reduce algorithm performance so it
#    is not recommended.  By default this step is not taken.  If you want to
#    use this step then set the adjLPE parameter "doMax" to TRUE.
#
# 2. The pi/2 correction of the variance is too large at small sample sizes.
#    The use of empirically estimated adjustments based on sample size is used
#    instead.  This increases the power of the test and generates a more
#    theoretically expected p value distribution.  The lpeAdj function
#    implements the use of the sample size based adjustments.  The user must
#    set the adjust1 and adjust2 variables to the sample size based values
#    found in the global variable ADJ.VALUES.  The empirical adjustment values
#    have been estimated for replicate sizes up to 10 only.
#
# The wrapper function runLPE does not use the first step and uses the second
# step by default.  This is the recommended manner in which LPEadj should be
# run.  One can use or not use step 1 and 2 independently though.  The
# functions myBaseOlig.error and lpeAdj can be run outside of the wrapper.

#
# Empirically based adjustments for 3-10 replicates, The correct adjustment can
# be found by using the replicate size as the index of the vector. The first
# two values are for sample sizes of 1 and 2.  In those cases LPE uses a
# different algorithm and no adjustment is made.
#
#ADJ.VALUES <- c(1, 1, 1.34585905516761 ,1.19363228146169 ,1.436849413109
#                  ,1.289652132873 ,1.47658053092781 ,1.34382984852146
#                  ,1.49972130857404, 1.3835405678718)
            

#
# Use LPE method to estimate fold change and variance. Note that there is a
# difference between LPE in bioconductor and LPE.test in Splus.  LPE sets
# all rank ordered variances up to the max to the max. The function
# baseOlig.error.step1 is recoded in this file so the resetting of the
# variance doesn't happen.  If doMax=T and doAdj=F then the original LPE
# algorithm is executed.
#
# dat - matrix of expression data of dimension n x 2k. There are n genes in
#      data with 2 conditions, and each condition has k reps. First k columns
#      will be condition 1 and last k columns will be condition 2.
# labels - vector containing column labels (groups) of dat.
#          eg. labels=c(0,0,0,1,1,1) indicates two groups of three replicates.
# doMax - boolean: if T then all variances below the max variance in the
#                  ordered distribution of variances are set to the maximum
#                  variance.  The default is False.  This step has been shown
#                  to reduce algorithm performance so it is recommended to keep
#                  this variable set to False.
#
# doAdj - run LPE with using variance adjustment value based on number of
#         replicates (hardcoded in adjValues) rather than pi/2.
# q -  is the quantile width; q=0.01 corresponds to 100 quantiles
#          i.e. percentiles. Bins/quantiles have equal number of genes and
#          are split according to the average intensity A
#
# Returns a data frame as described in the lpe documentation.
#
lpeAdj <- function(dat, labels=NULL, doMax=FALSE, doAdj=TRUE, q=.01) {

    if(is.null(dat)) {
        print("ERROR: data is null object")
        stop() 
    }
  
    p <- dim(dat)[1]            # number of probesets
    cols <- dim(dat)[2]         # number of chips in experiment


    # get number of replicates in both groups
    if(!is.null(labels)) {
  
        labs <- unique(labels)
        if(length(labs) !=2) {
           cat("ERROR: There must be two groups in data. ", length(labs),
               "groups found.\n")
           stop()
        }

        groupIndex <- labels == labs[1]
        dat <- cbind(dat[,groupIndex], dat[,!groupIndex])

        rep1 <- sum(labels == labs[1])
        rep2 <- sum(labels == labs[2])

        # make sure labels matches the data matrix in size
        if(rep1 + rep2 != cols) {
            print("ERROR: labels parameter does not match size of data")
            stop()
        }
    }
    else { 
        print("ERROR: labels is NULL.")
        stop()
    }

    # empirically based adjustments for 3-10 replicates, The first two values
    # are for reps 1 and 2.  In those cases LPE uses a different algorithm
    # and no adjustment is made.  The correct adjustment value can be fetched
    # using the replicate number as an index.
    adjValues <- c(1, 1, 1.34585905516761 ,1.19363228146169 ,1.436849413109
                   ,1.289652132873 ,1.47658053092781 ,1.34382984852146
                   ,1.49972130857404, 1.3835405678718)

    # If there are more than 10 replicates and doAdj is true then still use
    # pi/2 as adjustment.
    if(doAdj){

        if((rep1 > 10) | (rep2 > 10)) {
           print("Warning: for sample sizes greater than 10 the default of 
                 pi/2 is used for the variance adjustment.")
        }
        
        adj1 <- ifelse(rep1 > 10, pi/2, adjValues[rep1])
        adj2 <- ifelse(rep2 > 10, pi/2, adjValues[rep2])
    }
       
    # calculate base line error distributions
    var1 <- adjBaseOlig.error(dat[,1:rep1], setMax1=doMax, q=q)
    var2 <- adjBaseOlig.error(dat[,(rep1+1):cols], setMax1=doMax, q=q)

    # run LPE
    if(doAdj) {
        cat("variance adjustment values used: group 1: ", adj1, " group 2",
            adj2, "\n")
        results <- calculateLpeAdj(dat[,1:rep1],dat[,(rep1+1):cols],var1,var2,
                                   probe.set.name=c(1:p), adjust1=adj1,
                                   adjust2=adj2)
    } else {
        results <- lpe(dat[,1:rep1],dat[,(rep1+1):cols],var1,var2,
                       probe.set.name=c(1:p))
    }

    return(results)
} 



#
# LPE function that applies LPE method but uses user defined variance
# adjustments (adjust1, adjust2) instead of pi/2
#
# x: Replicated data from first experimental condition (as matrix 
#          or data-frame)
# y: Replicated data from second experimental condition (as matrix
#           or data-frame)
#
# basevar.x: Baseline distribution of first condition obtained from 
#          function baseOlig.error
#
# basevar.y: Baseline distribution of second condition obtained from 
#          function baseOlig.error
#
# df: Degrees of freedom used in fitting smooth.spline to estimates
#          of var.M for bins in A
#
# array.type: Currently supports oligo arrays
# probe.set.name: Gene IDs. By default if they are not provided then
#          1,2,3,... is assigned as GeneID
#
# trim.percent: Percent of (A, var.M) estimates to trim from low  end of A
#
calculateLpeAdj <- function (x, y, basevar.x, basevar.y, df = 10,
                             array.type = "olig", 
                   probe.set.name = "OLIG.probe.name", trim.percent = 5,
                    adjust1=1.57, adjust2=1.57) 
{

    n1 <- ncol(x)
    n2 <- ncol(y)
    ngenes <- nrow(x)
    if (n1 < 2 | n2 < 2) {
        stop("No replicated arrays!")
    }
    if (n1 > 2 | n2 > 2) {
        var.x <- basevar.x[, 2]
        var.y <- basevar.y[, 2]
        median.x <- basevar.x[, 1]
        median.y <- basevar.y[, 1]
        median.diff <- median.x - median.y
        std.dev <- sqrt(adjust1 *(var.x/n1) + adjust2*(var.y/n2))
        z.stats <- median.diff/std.dev
        data.out <- data.frame(x = x, median.1 = median.x, std.dev.1 = sqrt(var.x), 
            y = y, median.2 = median.y, std.dev.2 = sqrt(var.y), 
            median.diff = median.diff, pooled.std.dev = std.dev, 
            z.stats = z.stats)
        row.names(data.out) <- probe.set.name
        return(data.out)
    }
    if (n1 == 2 & n2 == 2) {
        var.x <- basevar.x[, 2]
        var.y <- basevar.y[, 2]
        median.x <- basevar.x[, 1]
        median.y <- basevar.y[, 1]
        median.diff <- median.x - median.y
        std.dev <- sqrt((var.x/n1) + (var.y/n2))
        z.stats <- median.diff/std.dev
        pnorm.diff <- pnorm(median.diff, mean = 0, sd = std.dev)
        p.out <- 2 * apply(cbind(pnorm.diff, 1 - pnorm.diff), 
            1, min)
        sf.xi <- smooth.spline(basevar.x[, 1], basevar.x[, 2], 
            df = df)
        var.x0 <- fixbounds.predict.smooth.spline(sf.xi, median.x)$y
        sf.xi <- smooth.spline(basevar.y[, 1], basevar.y[, 2], 
            df = df)
        var.y0 <- fixbounds.predict.smooth.spline(sf.xi, median.y)$y
        flag <- matrix(".", ngenes, 2)
        p.val <- matrix(NA, ngenes, 2)
        x.stat <- abs(x[, 1] - x[, 2])/sqrt(2 * var.x0)
        p.val[, 1] <- 2 * (1 - pnorm(x.stat))
        flag[p.val[, 1] < 0.01, 1] <- "*"
        flag[p.val[, 1] < 0.005, 1] <- "**"
        flag[p.val[, 1] < 0.001, 1] <- "***"
        y.stat <- abs(y[, 1] - y[, 2])/sqrt(2 * var.y0)
        p.val[, 2] <- 2 * (1 - pnorm(y.stat))
        flag[p.val[, 2] < 0.01, 2] <- "*"
        flag[p.val[, 2] < 0.005, 2] <- "**"
        flag[p.val[, 2] < 0.001, 2] <- "***"
        data.out <- data.frame(x = x, median.1 = median.x, std.dev.1 = sqrt(var.x), 
            p.outlier.x = p.val[, 1], flag.outlier.x = flag[, 
                1], y = y, median.2 = median.y, std.dev.2 = sqrt(var.y), 
            p.outlier.y = p.val[, 2], flag.outlier.y = flag[, 
                2], median.diff = median.diff, pooled.std.dev = std.dev, 
            z.stats = z.stats)
        row.names(data.out) <- probe.set.name
        return(data.out)
    }
  }



#
# setMax - if T then all variances below the max variance in the ordered
#          distribution of variances are set to the maximum variance.
#
adjBaseOlig.error <-function (y, stats = median, q = 0.01, min.genes.int = 10,
                             div.factor = 1, setMax1=FALSE)  {
    baseline.step1 <- adjBaseOlig.error.step1(y, stats = stats, setMax=setMax1,
                                             q=q)
    baseline.step2 <- adjBaseOlig.error.step2(y, stats = stats, setMax=setMax1,
        baseline.step1, min.genes.int = min.genes.int, div.factor = div.factor)
    return(baseline.step2)
}



#
# This is a function of LPE which is redefined here.
#
# setMax - if T then all variances below the max variance in the ordered
#          distribution of variances are set to the maximum variance.
#
adjBaseOlig.error.step1 <- function (y, stats = median, setMax=FALSE, q = 0.01,
                                    df = 10) {   
    
    AM <- am.trans(y)
    A <- AM[, 1]
    M <- AM[, 2]
    median.y <- apply(y, 1, stats, na.rm=TRUE)
    quantile.A <- quantile(A, probs = seq(0, 1, q), na.rm = TRUE)
    quan.n <- length(quantile.A) - 1
    var.M <- rep(NA, length = quan.n)
    medianAs <- rep(NA, length = quan.n)
    if (sum(A == min(A, na.rm=TRUE), na.rm=TRUE) > (q * length(A))) {
        tmpA <- A[!(A == min(A, na.rm=TRUE))]
        quantile.A <- c(min(A, na.rm=TRUE), quantile(tmpA, probs = seq(q, 
            1, q), na.rm = TRUE))
    }
    for (i in 2:(quan.n + 1)) {
        # get number of variances for this quantiles As
        n.i <- length(!is.na(M[A > quantile.A[i - 1] & A <= quantile.A[i]]))
        # calculate scaling factor due to taking minus both ways
        mult.factor <- 0.5 * ((n.i - 0.5)/(n.i - 1))
        # variance of variances for this quantiles As
        var.M[i - 1] <- mult.factor * var(M[A > quantile.A[i - 
            1] & A <= quantile.A[i]], na.rm = TRUE)
        # get median of this quantiles As
        medianAs[i - 1] <- median(A[A > quantile.A[i - 1] & A <= 
            quantile.A[i]], na.rm = TRUE)
    }

    # DANGER DANGER DANGER
    # this line sets all variances below the maximum variance to the maximum.
    # This has been shown to be detrimental to algorithm performance so it is
    # recommended to keep setMax to the default position of False.
    if(setMax) {
       print("adjusting variance badly")
       var.M[1:which(var.M == max(var.M))] <- max(var.M)
    }

    base.var <- cbind(A = medianAs, var.M = var.M)
    sm.spline <- smooth.spline(base.var[, 1], base.var[, 2], df = df)
    min.Var <- min(base.var[, 2])
    var.genes <- fixbounds.predict.smooth.spline(sm.spline, median.y)$y
    var.genes[var.genes < min.Var] <- min.Var
    basevar.step1 <- cbind(A = median.y, var.M = var.genes)
    ord.median <- order(basevar.step1[, 1])
    var.genes.ord <- basevar.step1[ord.median, ]
    return(var.genes.ord)
}


adjBaseOlig.error.step2 <- function (y, baseOlig.error.step1.res, df = 10, 
                                    stats = median, setMax=FALSE,
                                    min.genes.int = 10, div.factor = 1) 
{
 
    AM <- am.trans(y)
    A <- AM[, 1]
    M <- AM[, 2]
    median.y <- apply(y, 1, stats, na.rm=TRUE)
    var.genes.ord <- baseOlig.error.step1.res
    genes.sub.int <- n.genes.adaptive.int(var.genes.ord,
                                          min.genes.int = min.genes.int, 
                                          div.factor = div.factor)
    j.start <- 1
    j.end <- 0
    var.M.adap <- rep(NA, length = length(genes.sub.int))
    medianAs.adap <- rep(NA, length = length(genes.sub.int))
    for (i in 2:(length(genes.sub.int) + 1)) {
        j.start <- j.end + 1
        j.end <- j.start + genes.sub.int[i - 1] - 1
        vect.temp<-(A> var.genes.ord[j.start, 1] & A <= var.genes.ord[j.end,1])
        n.i <- length(!is.na(M[vect.temp]))
        mult.factor <- 0.5 * ((n.i - 0.5)/(n.i - 1))
        var.M.adap[i - 1] <- mult.factor * var(M[vect.temp], 
            na.rm = TRUE)
        medianAs.adap[i - 1] <- median(A[vect.temp], na.rm = TRUE)
    }

    # This line will set all low intensity variances to the maximum variance if
    # setMax equals true.  It is recommended to not enable this practice unless
    # one is sure that it will affect only a small number of genes.  In certain
    # cases this practice can set the variances of a large number of genes to
    # the maximum which has a detrimental effect on detecting effects.
    if(setMax) {
       var.M.adap[1:which(var.M.adap == max(var.M.adap))] <- max(var.M.adap)
    }
    
    base.var.adap <- cbind(A.adap = medianAs.adap, var.M.adap = var.M.adap)
    sm.spline.adap <- smooth.spline(base.var.adap[, 1], base.var.adap[, 
        2], df = df)
    min.Var <- min(base.var.adap[, 2])
    var.genes.adap <- fixbounds.predict.smooth.spline(sm.spline.adap, 
        median.y)$y
    var.genes.adap[var.genes.adap < min.Var] <- min.Var
    basevar.all.adap <- cbind(A = median.y, var.M = var.genes.adap)
    return(basevar.all.adap)
}



