### R code from vignette source 'LPEadj.Rnw'

###################################################
### code chunk number 1: LPEadj.Rnw:259-266
###################################################
 # Loading the library and null dataset (two groups with three
 # replicates each).  
 library(LPEadj)
 dat <- matrix(rnorm(6000), ncol=6)
 
 # Applying LPE
 lpe.result <- lpeAdj(dat, labels=c(0,0,0,1,1,1), doMax=FALSE, doAdj=TRUE)


###################################################
### code chunk number 2: LPEadj.Rnw:275-294
###################################################
# Loading the library and null dataset (two groups with three
 # replicates each)
 library(LPEadj)
 dat <- matrix(rnorm(6000), ncol=6)

 ADJ.VALUES <- c(1, 1, 1.34585905516761 ,1.19363228146169 ,1.436849413109
                  ,1.289652132873 ,1.47658053092781 ,1.34382984852146
                  ,1.49972130857404, 1.3835405678718)
            
 # calculate base line error distributions
 var1 <- adjBaseOlig.error(dat[,1:3], setMax1=FALSE, q=.05)
 var2 <- adjBaseOlig.error(dat[,4:6], setMax1=FALSE, q=.05)

 # The correct variance adjustments can be fetched using the replicate 
 # number for each group as in index for the ADJ.VALUES vector.
 # eg: ADJ.VALUES[n] if there are n replicates in a group 
 results <- calculateLpeAdj(dat[,1:3],dat[,4:6],var1,var2,
                  probe.set.name=c(1:1000), adjust1=ADJ.VALUES[3],
                  adjust2=ADJ.VALUES[3])


