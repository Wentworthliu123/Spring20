#### Appplied Multivariate Analysis
#### A collection of the R codes used in Bus41912 over the years.
#### Putting together on June 14, 2014 by Ruey S. Tsay
####
"Behrens" <- function(x1,x2){
# The x1 and x2 are two data matrices with xi for population i.
# Test the equal mean vectors.
# Written by Ruey S. Tsay on April 17, 2008
if(!is.matrix(x1))x1=as.matrix(x1)
if(!is.matrix(x2))x2=as.matrix(x2)
n1 = dim(x1)[1]
n2 = dim(x2)[1]
p1 = dim(x1)[2]
p2 = dim(x2)[2]
if (p1 == p2){
x1bar=matrix(colMeans(x1),p1,1)
x2bar=matrix(colMeans(x2),p2,1)
dev = x1bar-x2bar
S1 = cov(x1)
S2 = cov(x2)
S1n = (1/n1)*S1
S2n = (1/n2)*S2
Sp =  S1n+S2n
Si=solve(Sp)
T2 = t(dev)%*%Si%*%dev
S1s= S1n%*%Si
S2s = S2n%*%Si
SS1s = S1s%*%S1s
SS2s = S2s%*%S2s
d1 = (sum(diag(SS1s))+(sum(diag(S1s)))^2)/n1 + (sum(diag(SS2s))+(sum(diag(S2s)))^2)/n2
v = (p1+p1^2)/d1
cat("Estimate of v: ",v,"\n")
deg = v-p1+1
tt=T2*deg/(v*p1)
pvalue=1-pf(tt,p1,deg)
Behrens=matrix(c(T2,pvalue),1,2)
colnames(Behrens) <- c("T2-stat","p.value")
}
cat("Test result:","\n")
print(Behrens,digit=4)
}

####
"BoxM" <- function(x,nv){
# The x is the data vector with the first n1 rows belonging to population 1
#  the (n1+1):(n1+n2) rows belonging to population 2, etc.
# nv = (n1,n2,...,ng)'
# The number of groups is the length of nv-vector.
# Box's M-test for equal covariance matrics
# Written by Ruey S. Tsay on April 18, 2008
Box.M=NULL
g=length(nv)
p=dim(x)[2]
S=array(0,dim=c(p,p,g))
Sp=matrix(0,p,p)
n=sum(nv)
deg2=n-g
M = 0
# tmp1 is the sum[(n_i-1)*ln(det(S_i))
# u1 is the sum[1/(n_i-1)]
tmp1=0
u1 = 0
idx=0
for (i in 1:g){
da=x[(idx+1):(idx+nv[i]),]
smtx=cov(da)
S[,,i]=smtx
Sp=(nv[i]-1)*smtx+Sp
tmp1=(nv[i]-1)*log(det(smtx))+tmp1
u1 = u1 + 1.0/(nv[i]-1)
print("determinant")
print(det(smtx))
idx=idx+nv[i]
}
Sp=Sp/deg2
M=deg2*log(det(Sp))-tmp1
u = (u1-(1.0/deg2))*(2*p^2+3*p-1)/(6*(p+1)*(g-1))
C = (1-u)*M
nu=p*(p+1)*(g-1)/2
pvalue=1-pchisq(C,nu)
Box.M=cbind(Box.M,c(C,pvalue))
row.names(Box.M)=c("Box.M-C","p.value")
cat("Test result:","\n")
print(Box.M)
BoxM <-list(Box.M=M, Test.Stat=C,p.value=pvalue)
}

####
"BoxCox" <- function(da,interval=c(-2,2)){
if(!is.matrix(da))da=as.matrix(da)
n = dim(da)[1]; k=dim(da)[2]
#### objective function 
 lklam <- function(lam,x){
   x = abs(x)
   if(min(x) < 0.0001){x = x+0.2}
   if(abs(lam) < 0.00001){
    y = log(x)
    }
    else{
     y = (x^lam-1)/lam
     }
   n1 = length(x)
   dev=scale(y,center=T,scale=F)
   s = sum(dev*dev)/n1
   nlk = -0.5*n1*log(s)+(lam-1)*sum(log(x))
   lklam = nlk
  }
 lambda = seq(min(interval),max(interval),0.05)
 est=NULL
 for (j in 1:k){
  fv=NULL; x=da[,j]
  nT=length(lambda)
  for (i in 1:nT){
   fv=c(fv,lklam(lambda[i],x))
   }
  idx = which.max(fv)
  est=c(est,lambda[idx])
  }
 cat("Estimated transformation: ",est,"\n")
}

####
"classify56" <- function(da,size,eqP=T,pr=c(0),newdata=NULL){
### classification using equation (11-56) of Johnson and Wickern.
###
# da: data matrix. The data are arranged according to populations.
# size: a vector of sample size of each populations so that 
#    the length g [g=length(size)] is the number of populations.
# eqP: switch for equal probabilities
# eqV: switch for equal covariance matrices.
## pr: prior probability
# Assume equal costs.
# Assume normality.
#
 if(!is.matrix(da))da <- as.matrix(da)
 nr = dim(da)[1]
 nc = dim(da)[2]
 nrnew <- 0; newclass <- NULL
 if(length(newdata) > 0){
  if(!is.matrix(newdata)) newdata <- as.matrix(newdata)
  nrnew=dim(newdata)[1]; ncnew=dim(newdata)[2]
  if(ncnew < nc){
    cat("newdata is not in the proper format","\n")
    return
  }
 }
if(is.null(newdata)){newdata=da; nrnew=nr; ncnew=nc}
##
 g=length(size)
# compute the sample mean and covariance matrix of each populations
 cm=matrix(0,nc,g)
 dm=g*nc
# S: stores the covariance matrices (one on top of the other)
## Sinv: stores the inverse covariance matrices
 S=matrix(0,dm,nc)
 Sinv=matrix(0,dm,nc)
 ist=0
 for (i in 1:g){
 x = da[(ist+1):(ist+size[i]),]
 if(nc > 1){
 cm1 = apply(x,2,mean)}
 else{
  cm1=c(mean(x))
  }
 cm[,i]=cm1
 smtx=var(x)
 jst =(i-1)*nc
 S[(jst+1):(jst+nc),]=smtx
 Si=solve(smtx)
 Sinv[(jst+1):(jst+nc),]=Si
 cat("Population: ",i,"\n")
 print("mean vector:")
 print(cm1,digits=3)
 print("Covariance matrix:")
 print(smtx,digits=4)
# print("Inverse of cov-mtx:")
# print(Si,digits=4)
 ist=ist+size[i]
}
##
if(eqP){
pr=rep(1,g)/g
}
##
 print("Assume equal covariance matrices and costs")
 Sp=matrix(0,nc,nc)
 for (i in 1:g){
 jdx=(i-1)*nc
 smtx=S[(jdx+1):(jdx+nc),]
 Sp=Sp+(size[i]-1)*smtx
}
 Sp = Sp/(nr-g)
 print("Sp")
 print(Sp,digits=4)
 Spinv=solve(Sp)
 print("Sp-inv")
 print(Spinv,digits=4)
##
### Classification
 for (i in 1:nrnew){
  xo=as.numeric(newdata[i,])
  xo <- matrix(xo,nc,1)
  for (k in 1:g){
##### compute the constant terms
   dki = rep(0,g)
   ckk=cm[,k]
   for (ii in 1:g){
    dif=matrix(c(ckk-cm[,ii]),1,nc); s=c(ckk+cm[,ii])
    dki[ii] = -0.5*dif%*%Spinv%*%matrix(s,nc,1)+dif%*%Spinv%*%xo
    dki[ii] = dki[ii]-log(pr[ii]/pr[k])
    }
    dki=dki[-k]
    if(min(dki) > 0){newclass=c(newclass,k)}
   }
 }
 cat("New classification: ","\n")
 print(newclass)
return <- list(newclass=newclass)
}

#####
"Cmeans" <- function(da,size,eqV=T,alpha=0.05){
# The data matrix is "da".
# size: the sample size for each population.
# eqV: equal covariance matrices.
# Written by Ruey S. Tsay for Bus 41912 class.
#
if(!is.matrix(da))da=as.matrix(da)
 nc=ncol(da)
 nr=nrow(da)
 n1=size[1]
 n2=size[2]
 mu1 = apply(da[1:n1,],2,mean)
 S1=cov(da[1:n1,])
 print("Population 1:")
cat("Mean-vector:","\n")
 print(mu1,digits=3)
cat("Covariance matrix","\n")
 print(S1,digits=3)
 mu2 = apply(da[(n1+1):nr,],2,mean)
 S2=cov(da[(n1+1):nr,])
 print("Population 2:")
cat("Mean-vector: ","\n")
 print(mu2,digits=3)
cat("Covariance matrix: ","\n")
 print(S2,digits=3)
 if(eqV){
 Sp = ((n1-1)/(nr-2))*S1+((n2-1)/(nr-2))*S2
 print("Pooled covariance matrix:")
 print(Sp,digits=3)
 Spm = ((1/n1) + (1/n2))*Sp
 Spmv = solve(Spm)
 }
else {
 Spm = S1*(1/n1)+S2*(1/n2)
 Spmv = solve(Spm)
 }
 dif=matrix((mu1-mu2),nc,1)
cat("differnces in means: ","\n")
 print(dif,digits=3)
 T2 = t(dif)%*%Spmv%*%dif
 d2=nr-nc-1
 csq=(nr-2)*nc/(nr-nc-1)
 cri=T2/csq
 p=1-pf(cri,nc,d2)
 print("Hotelling T2, approx-F, & its p-value")
 tmp=c(T2,cri,p)
 print(tmp)
# C.I. for differences in means
 print("Simultaneous Tsq. C.I. for difference in means")
 pr=1-alpha
 c = sqrt(csq*qf(pr,nc,d2))
 for (i in 1:nc){
 tmp=sqrt(Spm[i,i])*c
 ci = c(dif[i]-tmp,dif[i]+tmp)
 print(ci,digits=3)
}
# Bonferroni simultanrous C.I. 
 print("Simultaneous Bonferroni C.I. for difference in means")
 pr=1 - alpha/(2*nc)
 cri=qt(pr,d2)
 for (i in 1:nc){
 tmp=cri*sqrt(Spm[i,i])
 ci = c(dif[i]-tmp,dif[i]+tmp)
 print(ci,digits=3)
}
# Print out the most responsible linear combination
 av=Spmv%*%dif
 print("Critical linear combination: ")
 print(av,digits=3)

}

### Given the data
"confreg" <- function(da,alpha=0.05,length=FALSE){
# The data matrix is "da".
# Written by Ruey S. Tsay for Bus 41912 class.
#
if(!is.matrix(da))da=as.matrix(da)
 nr = nrow(da)
 nc = ncol(da)
 cm=matrix(colMeans(da),1,nc)
 s=cov(da)
 simucr=matrix(0,nc,2)
 dg2=nr-nc
 cr=qf((1-alpha),nc,dg2)
 cr1=sqrt(nc*(nr-1)*cr/(nr-nc))
 se=sqrt(diag(s))/sqrt(nr)
 for (i in 1:nc){
 simucr[i,1]=cm[i]-cr1*se[i]
 simucr[i,2]=cm[i]+cr1*se[i]
}
 print("C.R. based on T^2")
 print(simucr)
 
 indvcr=matrix(0,nc,2)
 q=1-(alpha/2)
 cr=qt(q,(nr-1))
 for (i in 1:nc){
 indvcr[i,1]=cm[i]-cr*se[i]
 indvcr[i,2]=cm[i]+cr*se[i]
}
 print("CR based on individual t")
 print(indvcr)

 bonfcr=matrix(0,nc,2)
 q=1-(alpha/(2*nc))
 cr=qt(q,(nr-1))
 for (i in 1:nc){
 bonfcr[i,1]=cm[i]-cr*se[i]
 bonfcr[i,2]=cm[i]+cr*se[i]
}
 print("CR based on Bonferroni")
 print(bonfcr)

 asymcr=matrix(0,nc,2)
 cr=sqrt(qchisq((1-alpha),nc))
 for (i in 1:nc) {
 asymcr[i,1]=cm[i]-cr*se[i]
 asymcr[i,2]=cm[i]+cr*se[i]
}
 print("Asymp. simu. CR")
 print(asymcr)

 if(length){
  print("Lengths of confidence intervals:")
  leng=matrix(0,nc,4)
  leng[,1]=simucr[,2]-simucr[,1]
  leng[,2]=indvcr[,2]-indvcr[,1]
  leng[,3]=bonfcr[,2]-bonfcr[,1]
  leng[,4]=asymcr[,2]-asymcr[,1]
  colnames(leng) <- c("T^2","ind-t","Bonf","Asym")
  print(leng,digits=3)
  }
}

####
"contrast" <- function(da,cmtx,alpha=0.05){
# The data matrix is "da".
# Written by Ruey S. Tsay for Bus 41912 class.
#
if(!is.matrix(da))da=as.matrix(da)
 nr = nrow(da)
 nc = ncol(da)
 ave=matrix(colMeans(da),1,nc)
 S=cov(da)
 cm = cmtx%*%t(ave)
 co=cmtx%*%S%*%t(cmtx)
 Si=solve(co)
 Tsq = nr*(t(cm)%*%Si%*%cm)
 qm1=dim(cmtx)[1]
 tmp = Tsq*(nr-qm1)/((nr-1)*qm1)
 deg2=nr-qm1
 pv=1-pf(tmp,qm1,deg2)
 print("Hotelling T^2 statistics & p-value")
 print(c(Tsq,pv))
 print("Simultaneous C.I. for each contrast")
 cri = sqrt((nr-1)*qm1*qf(1-alpha,qm1,deg2)/deg2)
 simucr=matrix(0,qm1,2)
 for (i in 1:qm1){
 c=cmtx[i,]
 me=t(c)%*%t(ave)
 s = t(c)%*%S%*%c
 se=sqrt(s)/sqrt(nr)
 simucr[i,1]=me-cri*se
 simucr[i,2]=me+cri*se
}
 print(simucr,digits=3)
}

#### Given sample mean and sample covariance ###
"confreg.s" <- function(sm,s,nr,alpha=0.05){
# sm is the sample mean vector
# s: sample covariance matrix
# nr: number of observations
# written by Ruey S. Tsay for Bus 41912 class.
if(!is.matrix(s))s=as.matrix(s)
 nc=ncol(s)
 cm=matrix(sm,1,nc)
 simucr=matrix(0,nc,2)
 dg2=nr-nc
 cr=qf((1-alpha),nc,dg2)
 cr1=sqrt(nc*(nr-1)*cr/(nr-nc))
 se=sqrt(diag(s))/sqrt(nr)
 for (i in 1:nc){
 simucr[i,1]=cm[i]-cr1*se[i]
 simucr[i,2]=cm[i]+cr1*se[i]
}
 print("C.R. based on T^2")
 print(simucr)

 indvcr=matrix(0,nc,2)
 q=1-(alpha/2)
 cr=qt(q,(nr-1))
 for (i in 1:nc){
 indvcr[i,1]=cm[i]-cr*se[i]
 indvcr[i,2]=cm[i]+cr*se[i]
}
 print("CR based on individual t")
 print(indvcr)

 bonfcr=matrix(0,nc,2)
 q=1-(alpha/(2*nc))
 cr=qt(q,(nr-1))
 for (i in 1:nc){
 bonfcr[i,1]=cm[i]-cr*se[i]
 bonfcr[i,2]=cm[i]+cr*se[i]
}
 print("CR based on Bonferroni")
 print(bonfcr)

 asymcr=matrix(0,nc,2)
 cr=sqrt(qchisq((1-alpha),nc))
 for (i in 1:nc) {
 asymcr[i,1]=cm[i]-cr*se[i]
 asymcr[i,2]=cm[i]+cr*se[i]
}
 print("Asymp. simu. CR")
 print(asymcr)

}

##############

######
"eigTest" <- function(Sigma,p,q,n){
### Perform eigenvalue test using asymptotic chi-squares
### Sigma: the (p+q)-by-(p+q) covariance or correlation matrix
### n: sample size
###
###
   if(!is.matrix(Sigma))Sigma=as.matrix(Sigma)
   Sig=(Sigma+t(Sigma))/2
   S11=Sig[1:p,1:p]
   S12=Sig[1:p,(p+1):(p+q)]
   S22=Sig[(p+1):(p+q),(p+1):(p+q)]
   if(q < p){
     tmp = S11; S11=S22; S22=tmp
     S12=t(S12)
     }
   S11inv=solve(S11)
   S22inv=solve(S22)
   S=S11inv%*%S12%*%S22inv%*%t(S12)
   m1=eigen(S)
   minpq=min(p,q)
   vec1 = m1$vectors
   vec2 = S22inv%*%t(S12)%*%vec1
   tmp=crossprod(vec2,vec2)
   se = sqrt(diag(tmp))
   se=diag(1/se)
   vec2=vec2%*%se
   tst = 0; result=NULL
   adj=-(n-1-0.5*(p+q+1))
   for (j in 1:minpq){
     k1 = minpq-j+1
     k=minpq-j
     tst = tst + adj*log(1-m1$values[k1])
     df = (p-k)*(q-k)
     pv=1-pchisq(tst,df)
     result=rbind(result,c(j,tst,df,pv))
    }
    colnames(result) <- c("N(zero-Eig)","Test-statistic","df","p-value")
    print(result)
    cat("First set of ordered canonical variates: ","\n")
    print(vec1,digits=3)
    cat("Second set of ordered canonical variates: ","\n")
    print(vec2,digits=3)
  eigTest <- list(values=m1$values, Xvectors=vec1, Yvectors=vec2)
 }
    
#####
"EMmiss" <- function(da,fix=NULL,iter=1){
## da: data
## fix: an indicator matrix of the same dimension as da.
## fix(i,j) = 0 indicates missing, = 1 means observed.
## iter: number of iterations
##
if(!is.matrix(da))da=as.matrix(da)
n = dim(da)[1]
k = dim(da)[2]
x=da
if(length(fix) < 1)fix=matrix(1,n,k)
if(!is.matrix(fix))fix=as.matrix(fix)
# Compute the sample mean and covariance matrix.
mu=rep(0,k)
for (i in 1:k){
mu[i]=crossprod(x[,i],fix[,i])/sum(fix[,i])
}
cat("Sample mean: ","\n")
print(mu)
## compute the sample covariance matrix
xx = x
for (i in 1:k){
idx=c(1:n)[fix[,i]==0]
x[idx,i]=mu[i]
xx[idx,i]=0
}
###print(x)
S=cov(x)*(n-1)/n
cat("Sample covariance matrix (MLE): ","\n")
print(S)
### Compute the observed sum of squares
T2=t(xx)%*%xx
### Perform EM iterations
for (it in 1:iter){
### locate missing values in each data point (row of x-matrix)
S2=T2
for (i in 1:n){
idx=c(1:k)[fix[i,]==0]
if(length(idx) > 0){
jdx=c(1:k)[-idx]
obs=x[i,jdx]
###cat("obs: ","\n")
###print(obs)
mu1=mu[idx]
mu2=mu[jdx]
S11=S[idx,idx]
S12=S[idx,jdx]
S21=S[jdx,idx]
S22=S[jdx,jdx]
S22i=solve(S22)
x1=mu1+S12%*%S22i%*%(obs-mu2)
x[i,idx]=x1
X11=S11-S12%*%S22i%*%S21+x1%*%t(x1)
X12=x1%*%t(obs)
##cat("X11 and X12","\n")
##print(X11)
##print(X12)
S2[idx,idx]=S2[idx,idx]+X11
S2[idx,jdx]=S2[idx,jdx]+X12
S2[jdx,idx]=t(S2[idx,jdx])
S2[jdx,jdx]=S2[jdx,jdx]
}
#
}
cat("Iteration: ",it,"\n")
##cat("Modified data matrix: ","\n")
##print(x)
## compute the mean and covariance matrix 
mu=colMeans(x)
cat("iterated Sample mean: ","\n")
print(mu)
## compute the sample covariance matrix
muv=matrix(mu,k,1)
####print(S2)
S=S2/n - muv%*%t(muv)
cat("Iterated sample covariance matrix (MLE): ","\n")
print(S)
}
## end of the program
}

#####
"growth" <- function(da,nv,tp,q){
# This program tests the q-order growth curve model.
# The da is the data matrix.
# nv is the vector of sample sizes of the groups.
# q is the polynomial order.
# tp: time vector
# Written by Ruey S. Tsay on April 24, 2008
growth=NULL
N = dim(da)[1]
p = dim(da)[2]
g=length(nv)
Sp=matrix(0,p,p)
idx=0
xbar=NULL
for (i in 1:g){
ni=nv[i]
x=da[(idx+1):(idx+ni),]
ave=apply(x,2,mean)
xbar=cbind(xbar,ave)
S=cov(x)
Sp=(ni-1)*S+Sp
idx=idx+ni
}
W=Sp
Sp=Sp/(N-g)
Si=solve(Sp)
beta=matrix(0,(q+1),g)
semtx=matrix(0,(q+1),g)
B=matrix(1,p,1)
for (j in 1:q){
B=cbind(B,tp^j)
}
idx=0
Wq=matrix(0,p,p)
k=(N-g)*(N-g-1)/((N-g-p+q)*(N-g-p+q+1))
for (i in 1:g){
ni=nv[i]
x=da[(idx+1):(idx+ni),]
S=t(B)%*%Si%*%B
Sinv=solve(S)
d1=xbar[,i]
xy=t(B)%*%Si%*%d1
bhat=Sinv%*%xy
beta[,i]=bhat
semtx[,i]=sqrt((k/ni)*diag(Sinv))
fit=t(B%*%bhat)
err=x-kronecker(fit,matrix(1,ni,1))
err=as.matrix(err)
Wq=t(err)%*%err+Wq
idx=idx+ni
}
Lambda=det(W)/det(Wq)
Tst=-(N-(p-q+g)/2)*log(Lambda)
pv=1-pchisq(Tst,(p-q-1)*g)
growth=cbind(growth,c(Tst,pv))
row.names(growth)=c("LR-stat","p.value")
print("Growth curve model")
print("Order: ")
print(q)
print("Beta-hat: ")
print(beta,digits=4)
print("Standard errors: ")
print(semtx,digits=4)
print("W")
print(W,digits=4)
print("Wq")
print(Wq,digits=4)
print("Lambda:")
print(Lambda)
print("Test result:")
print(growth,digits=4)
}

####
"Hotelling" <- function(da,mu=NULL){
# The data matrix is "da".
# The mean vector is mu. (a list of numbers).
#
if(!is.matrix(da))da=as.matrix(da)
Hotelling=NULL
nr = dim(da)[1]
nc = dim(da)[2]
if(is.null(mu))mu = matrix(rep(0,nc),nc,1)
cm=matrix(colMeans(da),nc,1)
S = cov(da)
si=solve(S)
mu0=matrix(mu,nc,1)
dev=cm-mu0
T2 = nr*(t(dev)%*%si%*%dev)
d2=nr-nc
tt=T2*d2/(nc*(nr-1))
pvalue=1-pf(tt,nc,d2)
Hotelling=cbind(Hotelling,c(T2,pvalue))
row.names(Hotelling)=c("Hoteliing-T2","p.value")
Hotelling
}

####
"mlrchk" <- function(z,y,constant=TRUE,pred=NULL,alpha=0.05){
# Perform multiple linear regression analysis, provide leverage and 
# influential information, and compute prediction intervals.
#
#  Created April 27, 2010. Ruey Tsay
#
if(!is.matrix(z))z=as.matrix(z)
n=dim(z)[1]
r=dim(z)[2]
y1=as.matrix(y,n,1)
# estimation
z1=z
if(constant)z1=cbind(rep(1,n),z)
r1=dim(z1)[2]
ztz=t(z1)%*%z1
zty=t(z1)%*%y
ztzinv=solve(ztz)
beta=ztzinv%*%zty
res=y-z1%*%beta
sig=sum(res^2)/(n-r1)
# print results
print("coefficient estimates:")
pmtx=NULL
for (i in 1:r1){
s1=sqrt(sig*ztzinv[i,i])
tra=beta[i]/s1
ii = i; if(constant)ii=ii-1
pmtx=rbind(pmtx,c(ii,beta[i],s1,tra))
}
colnames(pmtx) <- c("Regor","Estimate","StdErr","t-ratio")
print(pmtx,digit=3)
# Compute the hat-matrix
par(mfcol=c(3,1))
H = z1%*%ztzinv%*%t(z1)
hii=diag(H)
ImH=rep(1,n)-hii
stres=res/sqrt(ImH*sig)
plot(stres,type='h',xlab='index',ylab='st-resi')
abline(h=c(0))
title(main='Studentized residuals')
plot(hii,type='h',xlab='index',ylab='h(i,i)',ylim=c(0,1))
lines(1:n,rep(1/n,n),lty=2,col="blue")
title(main='Leverage plot')
# Compute the Cook's distance.
# starts with studentized residuals
CookD=(hii/ImH)*(stres^2)/r1
plot(CookD,type='h',xlab='index',ylab='Cook Dist')
title(main='Cook Distance')
# Compute the predictive interval
if(length(pred)>0){
z0=pred
if(constant)z0=c(1,z0)
z0=matrix(z0,r1,1)
fst=t(z0)%*%beta
tc=qt((1-alpha/2),(n-r1))
sd=t(z0)%*%ztzinv%*%z0*sig
sd1=sqrt(sd)
cint=c(fst-tc*sd1,fst+tc*sd1)
sd2=sqrt(sd+sig)
pint=c(fst-tc*sd2,fst+tc*sd2)
print("100(1-alpha)% confidence interval:")
print(cint,digit=3)
print("100(1-alpha)% predictive interval:")
print(pint,digit=3)
}
mlrchk <- list(Cook=CookD,Hii=hii,stres=stres)
}

####
"mmlr" <- function(y,z,constant=TRUE){
# This program performs multivariate multiple linear regression analysis.
# z: design matrix 
# constant: switch for the constant term of the regression model
# y: dependent variables
## Model is y = z%*%beta+error
if(!is.matrix(y)) y <- as.matrix(y)
m <- ncol(y)
##
z=as.matrix(z)
n=nrow(z)
nx=ncol(z)
zc=z
if (constant) zc=cbind(rep(1,n),z)
p=ncol(zc)
y=as.matrix(y)
ztz=t(zc)%*%zc
zty=t(zc)%*%y
ZtZinv=solve(ztz)
beta=ZtZinv%*%zty
Beta=beta
cat("Beta-Hat matrix: ","\n")
print(round(Beta,3))
res=y-zc%*%beta
sig=t(res)%*%res/(n-p)
cat("LS residual covariance matrix: ","\n")
print(round(sig,3))
co=kronecker(sig,ZtZinv)
cat("Individual LSE of the parameter","\n")
##print("  est   s.d.   t-ratio    prob")
par=beta[,1]
deg=n-p
p1=ncol(y)
if (p1 > 1){
for (i in 2:p1){
par=c(par,beta[,i])
}
}
iend=nrow(beta)*ncol(beta)
tmp=matrix(0,iend,4)
for (i in 1:iend){
sd=sqrt(co[i,i])
tt=par[i]/sd
pr=2*(1-pt(abs(tt),deg))
tmp[i,]=c(par[i],sd,tt,pr)
}
colnames(tmp) <- c("Estimate","stand.Err","t-ratio","p-value")
#print(tmp,digits=3)
print(round(tmp,3))
# testing
sigfull=sig*(n-p)/n
det1=det(sigfull)
cat("===================","\n")
cat("Test for overall mmlr: ","\n")
C0 <- cov(y)*(n-1)/n
det0 <- det(C0)
tst <- -(n-p-0.5*(m-nx+1))*(log(det1)-log(det0))
df <- m*nx
pv <- 1-pchisq(tst,df)
cat("Test statistic, df,  and p-value: ",c(tst,df, pv),"\n")
#
if (nx > 1){
ztmp=z
cat("===================","\n")
print("Testing individual regressor")
for (j in 1:nx){
zc <- ztmp[,-j]
if(constant) zc=cbind(rep(1,n),zc)
ztz=t(zc)%*%zc
zty=t(zc)%*%y
ztzinv=solve(ztz)
beta=ztzinv%*%zty
res1=y-zc%*%beta
sig1=t(res1)%*%res1/n
det2=det(sig1)
tst=log(det1)-log(det2)
tst=-(n-p-0.5*ncol(y))*tst
deg=ncol(y)
pr=1-pchisq(tst,deg)
tmp=matrix(c(j,tst,pr),1,3)
colnames(tmp) <- c("regressor", "test-stat", "p-value")
print(round(tmp,4))
}

}
mmlr <- list(beta=Beta,residuals=res,sigma=sig,ZtZinv=ZtZinv,y=y,z=z,intercept=constant)
}

##############
"mmlrTest" <- function(z1,z2,y,constant=TRUE){
# This program performs likelihood ratio test for multivariate multiple linear regression models. 
# z1: design matrix of the reduced model
# z2: the regressors to be tested. That is, the full model contains regressors[z1,z2]
# constant: switch for the constant term of the regression model
# y: dependent variables
if(!is.matrix(z1))z1=as.matrix(z1)
if(!is.matrix(z2))z2=as.matrix(z2)
if(!is.matrix(y))y=as.matrix(y)
# sample size
n=nrow(z1)
q=ncol(z1)
# estimate the full model
z=cbind(z1,z2)
m=ncol(y)
r=ncol(z)
if(constant) z <- cbind(rep(1,n),z)
ztz=t(z)%*%z
zty=t(z)%*%y
ztzinv=solve(ztz)
beta=ztzinv%*%zty
Betafull=beta
res=y-z%*%beta
sig=t(res)%*%res/n 
print("Beta-Hat of Full model:")
print(beta,digits=3)
# Estimate the reduced model
z=z1
if(constant)z=cbind(rep(1,n),z)
ztz=t(z)%*%z
zty=t(z)%*%y
ztzinv=solve(ztz)
beta=ztzinv%*%zty
BetaR=beta
res=y-z%*%beta
sig1=t(res)%*%res/n 
print("Beta-Hat of Reduced model:")
print(beta,digits=3)
d1=det(sig)
d2=det(sig1)
Test=-(n-r-1-0.5*(m-r+q+1))*(log(d1)-log(d2))
deg=m*(r-q)
pvalue=1-pchisq(Test,deg)
print("Likelihood ratio test: beta2 = 0 vs beta_2 .ne. 0")
print("Test & p-value:")
print(c(Test,pvalue))
print("degrees of freedom:")
print(deg)
# Other tests
E=n*sig
H=n*(sig1-sig)
Wilk=det(E)/det(E+H)
HEinv=solve(H+E)
Pillai=sum(diag(H%*%HEinv))
print("Wilk & Pillai statistics:")
print(c(Wilk,Pillai))

mmlrTest <- list(betaF=Betafull,betaR=BetaR,sigmaF=sig,sigmaR=sig1,LR=Test,Degrees=deg)
}

### mmlr confidence and prediction intervals
"mmlrInt" <- function(model,newx,alpha=0.05){
### newx should include 1 if the model has a constant term.
###
y <- model$y; z <- model$z; constant <- model$intercept
ZtZinv <- model$ZtZinv; beta <- model$beta; sigma=model$sigma
if(!is.matrix(newx))newx <- matrix(newx,1,length(newx))
r <- ncol(z)
n <- nrow(y)
if((ncol(newx)==r) && constant){
 newx <- cbind(rep(1,nrow(newx)),newx)
 }
yhat <- newx%*%beta
##
m <- ncol(y)
p <- ncol(newx)
npt <- nrow(newx)
fv <- qf(1-alpha,m,(n-r-m))
crit <- m*(n-p)/(n-r-m)
crit <- sqrt(crit*fv)
for (i in 1:npt){
cat("at predictors: ",newx[1,],"\n")
cat("Point prediction: ","\n")
print(round(yhat[i,],3))
wk <- matrix(c(newx[i,]),1,p)
tmp <- wk%*%ZtZinv%*%t(wk)
tmp1 <- tmp+1
confi <- NULL
predi <- NULL
for (j in 1:m){
  d1 <- tmp*sigma[j,j]*n/(n-p)
  d2 <- tmp1*sigma[j,j]*n/(n-p)
  lcl <- yhat[i,j]-crit*sqrt(d1)
  ucl <- yhat[i,j]+crit*sqrt(d1)
  lpi <- yhat[i,j]-crit*sqrt(d2)
  upi <- yhat[i,j]+crit*sqrt(d2)
  confi <- rbind(confi,c(lcl,ucl))
  predi <- rbind(predi,c(lpi,upi))
 }
cat("Simultaneous C.I. with prob",(1-alpha),"\n")
print(round(confi,4))
cat("Simultaneous P.I. with prob",(1-alpha),"\n")
print(round(predi,4))
}

}



#####
"profile" <- function(x1,x2){
# The x1 and x2 are two data matrices with xi for population i.
# This program performs profile analysis (Parallel, Coincident, Level).
# Written by Ruey S. Tsay on April 23, 2008
profile=NULL
n1 = nrow(x1)
n2 = nrow(x2)
p1 = ncol(x1)
p2 = ncol(x2)
if (p1 == p2){
x1bar=matrix(colMeans(x1),p1,1)
x2bar=matrix(colMeans(x2),p2,1)
dev = x1bar-x2bar
S1 = cov(x1)
S2 = cov(x2)
Sp =  ((n1-1)/(n1+n2-2))*S1+((n2-1)/(n1+n2-2))*S2
cmtx=matrix(0,(p1-1),p1)
for (i in 1:(p1-1)){
cmtx[i,i]=-1
cmtx[i,(i+1)]=1
}
print("Are the profiles parallel?")
d1=cmtx%*%dev
deg1=p1-1
deg2=n1+n2-p1
S=((1/n1)+(1/n2))*cmtx%*%Sp%*%t(cmtx)
Si=solve(S)
Tsq=t(d1)%*%Si%*%d1
c=Tsq*deg2/((n1+n2-2)*deg1)
pv=1-pf(c,deg1,deg2)
profile=cbind(profile,c(Tsq,pv))
row.names(profile)=c("Test-T2","p.value")
print(profile,digits=4)
print("Are the profiles coincident?")
one=matrix(1,p1,1)
d1=sum(dev)
s=t(one)%*%Sp%*%one
s=s*((1/n1)+(1/n2))
T2=d1^2/s
pv2=1-pf(T2,1,(n1+n2-2))
profile=cbind(profile,c(T2,pv2))
print(profile,digits=4)
print("Are the profiles level?")
x=t(cbind(t(x1),t(x2)))
xbar=matrix(colMeans(x),p1,1)
S=cov(x)
d1=cmtx%*%xbar
CSC=cmtx%*%S%*%t(cmtx)
Si=solve(CSC)
Tst=(n1+n2)*t(d1)%*%Si%*%d1
c3=Tst*(n1+n2-p1+1)/((n1+n2-1)*(p1-1))
pv3=1-pf(c3,deg1,(n1+n2-p1+1))
profile=cbind(profile,c(Tst,pv3))
print(profile,digits=4)
}
cat("Summary: ","\n")
profile
}

############
"profileM" <- function(X,loc=1){
### Profile analysis for more than two pupolations.
### X: data matrix, including group indicators
### loc: column number for the group indicator
if(!is.matrix(X))X=as.matrix(X)
X1=as.matrix(X[,-loc])
p = dim(X1)[2]  ## number of variables
if(p > 1){
  C = matrix(0,p-1,p)
  for (i in 1:(p-1)){
   C[i,i] = -1; C[i,i+1]=1
  }
  Y = X1%*%t(C)
  fac1 = factor(X[,loc])
  m1 = manova(Y~fac1)
  cat("Are the profiles parallel?","\n")
  print(summary(m1,test="Wilks"))
   cat("===","\n")
   cat("Are the profile coincident?","\n")
   J=matrix(1,p,1)
   Y1=X1%*%J
   m2=aov(Y1~fac1)
   print(summary(m2))
###
  cat("====","\n")
  cat("Are the profiles level? ","\n")
  Hotelling(Y,rep(0,p-1))
 }
}

#####
"qqbeta" <- function(da){
# The data matrix is "da".
 if(!is.matrix(da))da=as.matrix(da)
 nr = dim(da)[1]; nc = dim(da)[2]
 dev=scale(da,center=T,scale=F)
 dev=as.matrix(dev)      
 s=cov(da)            
 si=solve(s)
 d2=sort(diag(dev%*%si%*%t(dev))) 
 d2=nr*d2/((nr-1)^2)
 a <- nc/2; b <- (nr-nc-1)/2
 alpha <- (nc-2)/(2*nc)
 beta <- (nr-nc-3)/(2*(nr-nc-1))
 mn = min(a,b,alpha,beta)
 if(mn > 0){
  prob=(c(1:nr)-alpha)/(nr-alpha-beta+1)
  q1=qbeta(prob,a,b)   
  plot(q1,d2,xlab='Quantile of beta-dist',ylab='d^2')

  fit = lsfit(q1,d2)
  fitted = d2-fit$residuals
  lines(q1,fitted)
  rq=cor(q1,d2)
  cat("correlation coefficient:",rq,"\n")
 }
 else{
  cat("Insufficient sample size")
  }
}

#####
"qqchi2" <- function(da){
# The data matrix is "da".
 if(!is.matrix(da))da=as.matrix(da)
 nr = dim(da)[1]
 nc = dim(da)[2]
 dev=scale(da,center=T,scale=F)
 dev=as.matrix(dev)      
 s=cov(da)            
 si=solve(s)
 d2=sort(diag(dev%*%si%*%t(dev))) 
 prob=(c(1:nr)-0.5)/nr
 q1=qchisq(prob,nc)   
 plot(q1,d2,xlab='Quantile of chi-square',ylab='d^2')

 fit = lsfit(q1,d2)
 fitted = d2-fit$residuals
 lines(q1,fitted)
 rq=cor(q1,d2)
 cat("correlation coefficient:",rq,"\n")
}

####
"sir" <- function(y,x,H){
# Performs the Sliced Inverse Regression analysis
# y: the data matrix, nT by 1, of the dependent variable
# x: the nT-by-p matrix of predictors.
# H: the number of slices.
#
z=as.matrix(x)
p=ncol(z)
nT=nrow(z)
m1=sort(y,index.return=T)
idx=m1$ix
if(H < 5)H=5
m=floor(nT/H)
nH=rep(m,H)
rem = nT-m*H
if(rem > 0){
for (i in 1:rem){
nH[i]=nH[i]+1
}
}
#print(nH)
xmean=apply(z,2,mean)
ones=matrix(rep(1,nT),nT,1)
xm=matrix(xmean,1,p)
Xdev=z-ones%*%xm
varX=cov(Xdev)
m1=eigen(varX)
P=m1$vectors
Dinv=diag(1/sqrt(m1$values))
Sighlf=P%*%Dinv%*%t(P)
X1=Xdev%*%Sighlf
# the sliced mean vectors of standardized predictors
EZ=matrix(0,H,p)
X1=X1[idx,]

iend=0
for (i in 1:H){
ist=iend+1
iend=iend+nH[i]
for (j in ist:iend){
EZ[i,]=EZ[i,]+X1[j,]
}
EZ[i,]=EZ[i,]/nH[i]
ist=iend
}

Msir=matrix(0,p,p)
for (i in 1:H){
tmp=matrix(EZ[i,],p,1)
Msir=Msir+(nH[i]/nT)*tmp%*%t(tmp)
}
m2=eigen(Msir)
eiv=m2$values
print(eiv,digits=3)
vec=Sighlf%*%m2$vectors
print(vec,digits=3)
tX=z%*%vec

# perform tests
result=matrix(0,p,3)
print("Testing:")
for (i in 1:p){
j=p-i+1
avl=mean(eiv[j:p])
tst = nT*i*avl
deg = i*(H-1)
pva=1-pchisq(tst,deg)
result[i,]=c(i,tst,pva)
}
print(result,digits=4)

list(values=m2$values,dir=vec,transX=tX)

}

######
"t2chart" <- function(da){
# The data matrix is "da".
# Written by Ruey S. Tsay on April 10, 2008.
if(!is.matrix(da))da=as.matrix(da)
 nr = nrow(da)
 nc = ncol(da)
 dev=scale(da,center=T,scale=F)
 dev=as.matrix(dev)      
 s=t(dev)%*%dev/(nr-1)             
 si=solve(s)
 t2=diag(dev%*%si%*%t(dev)) 
 ul1=qchisq(0.95,nc)
 ul2=qchisq(0.99,nc)
 yl=max(t2,ul1,ul2)
 plot(t2,type='l',ylim=c(0,yl))
 points(t2,pch='*')
 abline(h=c(ul1),lty=2)
 abline(h=c(ul2),lty=3)
 title(main='The limits are 95% and 99% quantiles')
 t2chart <- list(Tsq=t2)
}

####
"t2future" <- function(da,ini){
# The data matrix is "da".
# Program written by Ruey S. Tsay, April 10, 2008.
# Modified April 15, 2014.
if(!is.matrix(da))da=as.matrix(da)
 nr = nrow(da)
 nc = ncol(da)
 if (ini < nc) ini= nc+2
 npts= nr-ini
 T2=rep(0,npts)
 for (i in ini:nr){
#### Data used is from t=1 to t=i-1.
  dat=da[1:(i-1),]
  dev=scale(dat,center=T,scale=F)
  cm=matrix(apply(dat,2,mean),1,nc)
  dev=as.matrix(dev) 
   ii=i-1     
   s=t(dev)%*%dev/(ii-1)
   si=solve(s)
   xi=matrix(da[i,],1,nc)
   t2= (ii/i)*(xi-cm)%*%si%*%t(xi-cm)
   ul1=((ii-1)*nc)/(ii-nc)*qf(0.95,nc,(ii-nc))
   ul2=((ii-1)*nc)/(ii-nc)*qf(0.99,nc,(ii-nc))
   T2[(i-ini+1)]=t2
 }
 yl=max(T2,ul2)
 tdx=c(ini:nr)
 plot(tdx,T2,type='l',ylim=c(0,yl))
 points(tdx,T2,pch='*')
 abline(h=c(ul1),lty=2)
 abline(h=c(ul2),lty=3)
 title(main='Future Chart: the limits are 95% and 99% quantiles')
 t2future <- list(Tsq=T2)
}

####
"NormSimMean" <- function(n,k=2,mean=rep(0,k),sigma=diag(rep(1,k)),iter=1){
## Use simulation to see the sampling properties of sample mean.
## n: sample size (can be a sequence of sample sizes)
## k: dimension
## iter: number of iterations in simulation
##
require(mvtnorm)
k1 <- length(mean)
k2 <- ncol(sigma)
k <- min(k1,k2)
##
nsize <- length(n)
sm <- NULL
scov <- NULL
for (it in 1:iter){
 for (j in 1:nsize){
  nob <- n[j]
  r <- rmvnorm(nob,mean=mean[1:k],sigma=sigma[1:k,1:k])
  mu <- apply(r,2,mean)
  v1 <- cov(r)
  sm <- rbind(sm,mu)
  scov <- rbind(scov,c(v1))
  }
 }
#
NormSimMean <- list(sm=sm,scov=scov)
}

"MeanSim" <- function(mu1,mu2,sigma1=NULL,sigma2=NULL,n1,n2,df=0,iter=1){
## Simualtion to compare mean vectors of two populations
##
require(mvtnorm)
k1 <- length(mu1); k2 <- length(mu2)
k <- min(k1,k2)
if(is.null(sigma1))sigma1 <- diag(rep(1,k))
if(is.null(sigma2))sigma2 <- diag(rep(1,k))
sigma1 <- (sigma1+t(sigma1))/2
sigma2 <- (sigma2+t(sigma2))/2
del <- mu1-mu2
#
if(n1 < 1)n1 <- 2*k
if(n2 < 1)n2 <- 2*k
fr <- (n1-1)/(n1+n2-2)
dmean <- NULL; Sigma1 <- NULL; Sigma2 <- NULL; Fstat <- NULL
for (it in 1:iter){
  if(df <= 0){
   x1 <- rmvnorm(n1,mean=mu1,sigma=sigma1)
   x2 <- rmvnorm(n2,mean=mu2,sigma=sigma2)
   }else{
   x1 <- rmvt(n1,sigma=sigma1,df=df)
   x1 <- x1 + matrix(1,n1,1)%*%matrix(mu1,1,k)
   x2 <- rmvt(n2,sigma=sigma2,df=df)
   x2 <- x2 + matrix(1,n2,1)%*%matrix(mu2,1,k)
   }
#
 sm1 <- apply(x1,2,mean)
 sm2 <- apply(x2,2,mean)
 S1 <- cov(x1)
 S2 <- cov(x2)
 dm <- sm1-sm2 - del
 Sp <- S1*fr+S2*(1-fr)
 Sp <- Sp*((1/n1)+(1/n2))
 Spinv <- solve(Sp)
 T2 <- matrix(dm,1,k)%*%Spinv%*%matrix(dm,k,1)
 dmean <- rbind(dmean,dm)
 Sigma1 <- rbind(Sigma1,c(S1))
 Sigma2 <- rbind(Sigma2,c(S2))
 Fst <- T2*(n1+n2-k-1)/((n1+n2-2)*k)
 Fstat <- c(Fstat,Fst)
}

df2 <- n1+n2-k-1
Fdist <- rf(10000,df1=k,df2=df2)
cat("F-stat is ditrtibuted as F with dfs: ",c(k,(n1+n2-k-1)),"\n")
par(mfcol=c(1,1))
qqplot(Fstat,Fdist)

MeanSim <- list(dmean = dmean,sigma1=Sigma1, sigma2=Sigma2,Fstat=Fstat)

}

### LASSO simulations
"LASSO.sim" <- function(n,p,coef=rep(0,p),sigma2=1.0,df=0,X=NULL, sigmaX = 1.0){
## Perform LASSO simualtion.
## n: sample size
## p: dimension
## coef: true coefficient vectors
## sigma2: noise varinace
## df = 0: Gaussian noise; df .ne. 0, then use t-noise
## X: the design matrix. It will be generated from Gaussian if not given.
## sigmaX: variance of columns of X
## 
## output: Y: response
## X: design matrix
## epsilon: true noise
##
if(is.null(X)){
npt <- n*p
X <- matrix(rnorm(npt)*sigmaX,n,p)
}
if(!is.matrix(X)) X <- as.matrix(X)
if(df==0){
 eps <- rnorm(n)*sqrt(sigma2)
 }
 else{
 eps <- rt(n,df=df)*sqrt(sigma2)
 }

Y <- X%*%matrix(coef,p,1)+eps

LASSO.sim <- list(Y=Y,X=X,epsilon=eps,beta=coef)
}

###
"CVlasso" <- function(x,y,sizesample,s=c(0.1,0.25,0.5,0.75,0.9),iter=10){
## Compute MSE of out-of-sample prediction of LASSO
## x: design matrix 
## y: dependent variable
## sizesample: subsample size
n1 <- length(s)
MSE <- matrix(0,n1,iter)
X <- x; Y <- y
for (i in 1:iter){
 train <- sample(1:nrow(X),sizesample)
 lasso.m <- lars(x=X[train,],y=Y[train],type="lasso")
 for (j in 1:n1){
  MSE[j,i] <- mean((predict(lasso.m,X[-train,],s=s[j],mode="fraction")$fit-Y[-train])^2)
  }
 }
 mse <- apply(MSE,1,mean)
 mse <- rbind(s,mse)
 rownames(mse) <- c("fraction","mean(MSE)")
 print(mse)
 x1 <- paste("s = ",s)
 tmse <- t(MSE)
 colnames(tmse) <- x1
 boxplot(tmse,ylab='MSE')
}

###
"CVglmnet" <- function(x,y,sizesample,sidx=c(5,10,20,30,40,50,60,70,80,90,95),iter=20,family="gaussian",type="link",alpha=1){
## Compute MSE of out-of-sample prediction of LASSO
## x: design matrix 
## y: dependent variable
## sizesample: subsample size
n1 <- length(sidx)
X <- x; Y <- y
MSE <- matrix(0,n1,iter)
sidx <- sidx/100
for (i in 1:iter){
 train <- sample(1:nrow(X),sizesample)
 glmnet.fit <- glmnet(x=X[train,],y=Y[train],family=family,alpha=alpha)
 iidx <- floor(sidx*length(glmnet.fit$lambda))
 s <-glmnet.fit$lambda[sidx]
 NX <- X[-train,]
 for (j in 1:n1){
  pp <- predict(glmnet.fit,newx=NX,s=s[j],type=type)
  MSE[j,i] <- mean((pp-Y[-train])^2)
  }
 }
 mse <- apply(MSE,1,mean)
 mse <- rbind(s,mse)
 rownames(mse) <- c("lambda","mean(MSE)")
 print(mse)
 x1 <- paste("s = ",sidx)
 tmse <- t(MSE)
 colnames(tmse) <- x1
 boxplot(tmse,ylab='MSE',main="glmnet")
}

###
"discrim" <- function(da,size,eqP=T,eqV=T,pr=c(0),newdata=NULL){
# da: data matrix. The data are arranged according to populations.
# size: a vector of sample size of each populations so that 
#    the length g [g=length(size)] is the number of populations.
# eqP: switch for equal probabilities
# eqV: switch for equal covariance matrices.
# Assume equal costs.
# Assume normality.
#
 da=as.matrix(da)
 nr = dim(da)[1]
 nc = dim(da)[2]
 nrnew=0
 if(length(newdata) > 0){
  newclass=NULL
  if(!is.matrix(newdata)){
   newdata <- as.matrix(newdata)
  }
  nrnew=dim(newdata)[1]; ncnew=dim(newdata)[2]
  if(ncnew < nc){
    cat("newdata is not in the proper format","\n")
    return
  }
 }
##
 g=length(size)
# compute the sample mean and covariance matrix of each populations
 cm=matrix(0,nc,g)
 dm=g*nc
# S: stores the covariance matrices (one on top of the other)
 S=matrix(0,dm,nc)
 Sinv=matrix(0,dm,nc)
 ist=0
 for (i in 1:g){
 x = da[(ist+1):(ist+size[i]),]
 if(nc > 1){
 cm1 = apply(x,2,mean)}
  else{
  cm1=c(mean(x))
   }
 cm[,i]=cm1
 smtx=var(x)
 jst =(i-1)*nc
 S[(jst+1):(jst+nc),]=smtx
 Si=solve(smtx)
 Sinv[(jst+1):(jst+nc),]=Si
 cat("Population: ",i,"\n")
 print("mean vector:")
 print(cm1,digits=3)
 print("Covariance matrix:")
 print(smtx,digits=4)
 print("Inverse of cov-mtx:")
 print(Si,digits=4)
 ist=ist+size[i]
}
##
if(eqP){
pr=rep(1,g)/g
}else{
 pr <- size/sum(size)
 }
##
if(eqV){
 print("Assume equal covariance matrices and costs")
 Sp=matrix(0,nc,nc)
 for (i in 1:g){
 jdx=(i-1)*nc
 smtx=S[(jdx+1):(jdx+nc),]
 Sp=Sp+(size[i]-1)*smtx
}
 Sp = Sp/(nr-g)
 print("Sp: pooled covariance matrix")
 print(Sp,digits=4)
 Spinv=solve(Sp)
 print("Sp-inv")
 print(Spinv,digits=4)
 if(g==2){
  d=cm[,1]-cm[,2]
  a=t(d)%*%Spinv
  y1bar=a%*%matrix(cm[,1],nc,1)
  y2bar=a%*%matrix(cm[,2],nc,1)
  mhat=(y1bar+y2bar)/2
  print("Linear classification function:")
  cat('coefficients: ',a,"\n")
  cat('m-hat:',mhat,"\n")
 }
# compute ai-vector and bi value
 ai=matrix(0,g,nc)
 bi=rep(0,g)
 for (i in 1:g){
 cmi = as.matrix(cm[,i],nc,1)
 ai[i,]= t(cmi)%*%Spinv
 bi[i]=(-0.5)*ai[i,]%*%cmi+log(pr[i])
}
 cat("Discriminate functions: matrix: ","\n")
  print(ai,digits=3)
 cat("Discriminate fucntions: constant: ","\n")
  print(bi,digits=3)
##
 ist = 0
 tmc=matrix(0,g,g)
# icnt: counts of misclassification
 icnt = 0
 mtab=NULL
 for(i in 1:g){
 cls=rep(0,g)
 jst=ist+1
 jend=ist+size[i]
 for (j in jst:jend){
 dix=bi
 for (ii in 1:g){
 dix[ii]=dix[ii]+sum(ai[ii,]*da[j,])
}
# Find the maximum
 kx=which.max(dix)
#
 cls[kx]=cls[kx]+1
 if(abs(kx-i)>0){
 icnt=icnt+1
# compute the posterior probabilities
 xo1=matrix(da[j,],1,nc)
 print(xo1)
 print(c(j,icnt))
 tmp=xo1%*%Spinv
 tmp1=sum(xo1*tmp)
 tmp1 = (tmp1 + det(Sp))/2
 pro = exp(dix-tmp1)
 pro=pro/sum(pro)
 case=c(j,i,kx,pro)
 mtab=rbind(mtab,case)
}
}
tmc[i,]=cls
ist=ist+size[i]
}
######
##### new data classification
 if(nrnew > 0){
  for (j in 1:nrnew){
   dix = bi
   for (ii in 1:g){
    dix[ii]=dix[ii]+sum(ai[ii,]*newdata[j,])
   }
  newclass=c(newclass,which.max(dix))
  }
 }
}
# Unequal covariance matrices
else {
 print("Assume Un-equal covariance matrices, but equal costs")
# The next loop computes the dix components that do not depend on the obs.
 bi=log(pr)
 for (i in 1:g){
 jst=(i-1)*nc
 Si=S[(jst+1):(jst+nc),]
 bi[i]=bi[i]-0.5*log(det(Si))
} 
 cat("Discriminate fucntions: constant: ","\n")
  print(bi,digits=3)
 ist = 0
 tmc=matrix(0,g,g)
# icnt: counts of misclassification
 icnt = 0
# mtab: stores the misclassification cases
 mtab=NULL
 for(i in 1:g){
 cls=rep(0,g)
# identify the cases for population i.
 jst=ist+1
 jend=ist+size[i]
 for (j in jst:jend){
# Compute the dix measure
 xo=matrix(da[j,],1,nc)
 dix=bi
 for (ii in 1:g){
 cmi = matrix(cm[,ii],1,nc)
 kst=(ii-1)*nc
 Si=Sinv[(kst+1):(kst+nc),]
 dev=xo-cmi
 tmp=as.matrix(dev)%*%Si
 dix[ii]=dix[ii]-0.5*sum(dev*tmp)
}
# Find the maximum
 kx=which.max(dix)
#
 cls[kx]=cls[kx]+1
 if(abs(kx-i)>0){
 icnt=icnt+1
# compute the posterior probabilities
 pro = exp(dix)
 pro=pro/sum(pro)
 case=c(j,i,kx,pro)
 mtab=rbind(mtab,case)
}
}
tmc[i,]=cls
ist=ist+size[i]
}
#### new data classification, if needed.
 if(nrnew > 0){
  for (j in 1:nrnew){
  dix=bi
  for (ii in 1:g){
   cmi = matrix(cm[,ii],1,nc)
   kst=(ii-1)*nc
   Si=Sinv[(kst+1):(kst+nc),]
   xo=matrix(newdata[j,],1,nc)
   dev=xo-cmi
   tmp=as.matrix(dev)%*%Si
   dix[ii]=dix[ii]-0.5*sum(dev*tmp)
   }
   newclass=c(newclass,which.max(dix))
  }
 }
}
# Output
if(nrow(mtab)>0){
print("Misclassification cases")
print("Case From  To    Post-prob")
print(mtab,digits=3)
}
print("Classification table:")
print(tmc)
print("Error rate:")
aper=icnt/nr
print(aper)

###
if(nrnew > 0){
cat("New data classification","\n")
print(newclass)
}

return <- list(newclass=newclass,misclassTable=tmc)
}

### This program performs canonical correlation test of independence
"cancorTest" <- function(x,y){
if(!is.matrix(x))x <- as.matrix(x)
if(!is.matrix(y))y <- as.matrix(y)
n1 <- nrow(x); p <- ncol(x)
n2 <- nrow(y); q <- ncol(y)
n <- min(n1,n2)
m1 <- cancor(x[1:n,],y[1:n,])
cor <- rev(m1$cor)
k <- min(p,q)
tmp <- log(1-cor^2)
tst <- cumsum(tmp)
tst <- -(n-1-0.5*(p+q+1))*tst
Test <- NULL
for (i in 1:k){
ii <- k-i
df <- (p-ii)*(q-ii)
pv <- 1-pchisq(tst[i],df)
Test <- rbind(Test,c(tst[i],df,pv))
}
colnames(Test) <- c("Test-stat","df","p-value")
print(round(Test,3))
cancorTest <- list(Test=Test)
}