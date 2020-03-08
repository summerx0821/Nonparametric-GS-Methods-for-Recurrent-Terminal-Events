
library(MASS)
library(gdata)
#setwd
source("TayobMurray.gstest.function.R")
example <- read.csv("example_data.csv")


#SET UP PARAMETERS
s=c(12,24,36,48)
Tau=12
space=1.5
alpha=0.05
#Prepare data used in the function
data1 <- list(E=example[example$group==1,]$E,
              X=example[example$group==1,]$X,
              delta=example[example$group==1,]$delta,
              Z_star=as.matrix(example[example$group==1,6:48]))
data2 <- list(E=example[example$group==2,]$E,
              X=example[example$group==2,]$X,
              delta=example[example$group==2,]$delta,
              Z_star=as.matrix(example[example$group==2,6:48]))


#OUR METHOD: OVERALL TAU RESTRICTED MEAN TEST
newRMST.summary1=TayobMurray.gstest.mean.cov(data1,s,Tau,space)
newRMST.summary2=TayobMurray.gstest.mean.cov(data2,s,Tau,space)
#Calculate test statistics T
n1=newRMST.summary1$NumberEntered
n2=newRMST.summary2$NumberEntered
norm.const=sqrt(n1*n2/(n1+n2))
Tstat=norm.const*abs(newRMST.summary1$RmeanST-newRMST.summary2$RmeanST)
Cov.Tstat=norm.const%*%t(norm.const)*(newRMST.summary1$cov+newRMST.summary2$cov)


#COMPARE WITH CRITICAL VALUES
v=s/max(s)
cum.alpha=2-2*pnorm(qnorm(1-alpha/2)/sqrt(v))
alpha.each=c(cum.alpha[1],(cum.alpha[-1]-cum.alpha[-length(s)])/(1-cum.alpha[-length(s)]))
Z.T=mvrnorm(1000000,rep(0,length(s)),Cov.Tstat)

critical.Z.T=Find.critical.Z(Z.T,alpha.each)
Reject.T=ifelse(Tstat>=critical.Z.T,1,0)


## Results:
# > critical.Z.T
# [1] 7.863780 5.826057 4.626109 3.987550
# > Tstat
# [1] 5.124873 8.708382 8.803949 9.461295
# > Reject.T
# [1] 0 1 1 1

