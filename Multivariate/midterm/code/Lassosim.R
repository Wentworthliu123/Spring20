### Simulation of LASSO regression
"Lassosim" <- function(n,p,sigma=1.0,s0=6,iseed=121){
# n: sample size
# p: dimension
# sigma: residual standard error
# s0: number of nonzero coefficients
#
set.seed(iseed)
if(p < 1)p=10
if(s0 > p)s0=p
if(n < 10) n=30
beta <- rep(0,p)
idx <- sample(1:p,s0,replace=FALSE)
idx <- sort(idx)
b <- runif(max=10,min=-10,s0)
beta[idx] <- b
truecoef=rbind(idx,b)
#
x <- matrix(rnorm(n*p),n,p)
yhat <- x%*%matrix(beta,p,1)
y <- yhat+sigma*rnorm(n)


return <- list(y=y,x=x,beta=beta,truecoef=truecoef)
} 

#### ell-2 simple linear regression boosting
"sBoot" <- function(y,X,v=0.01,m=1000){
### Input: y--- dependent variable
###        X: matrix of predictors
###        v: boosting rate (0 < v <= 1)
###        m: boosting iterations
###
if(!is.matrix(X))X <- as.matrix(X)
k <- ncol(X)
n <- length(y)
n <- min(n,nrow(X))
y <- y[1:n]
X <- X[1:n,]
if(v > 1)v <- 1
if(v <= 0)v <- 0.01
### removing means
 y <- scale(y,center=TRUE,scale=FALSE)
 X <- scale(X,center=TRUE,scale=FALSE)
###
yhat <- rep(0,n)
beta <- rep(0,k)
y1 <- y
for (iter in 1:m){
 c1 <- cor(y1,X)
 idx <- which.max(abs(c1))
 m1 <- lm(y1~-1+X[,idx])
 yhat <- yhat+v*m1$fitted.values
 y1 <- y1-v*m1$fitted.values
 beta[idx] <- beta[idx]+v*m1$coefficients
}

resi <- y-yhat
idx <- c(1:k)[abs(beta) > 0.0000001]
icnt <- length(idx)
### 
cat("ell-2 boosting via simple linear regression: ","\n")
cat("Both Y and X are mean-adjusted.","\n")
cat("v and m: ",c(v,m),"\n")
cat("number of predictors selected and number of predictors: ",c(icnt,k),"\n")

sBoot <- list(beta=beta,residuals=resi,m=m,v=v,selection=idx,count=icnt,yhat=yhat)
}