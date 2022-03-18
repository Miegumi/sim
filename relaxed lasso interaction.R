#模拟数据
library(doParallel)
library(glmnet)
library(MASS)
library(dplyr)
library(stringr)
library(magrittr)
library(tidyr)

  set.seed(12345)
  registerDoParallel(cores=4)
  n = 50 #data size
  p = 10 #data dim
  q=5  
  r=0    
  r.square<-0.8  
  mu<- rep(0,p)  
  Sigma <- matrix(NA, ncol=p,nrow=p)  
  for(i in 1:p){ 
    for(j in 1:p){
      if(i!=j){Sigma[i,j]=r^abs(i-j)}  
      else{Sigma[i,j]=1}  
    }
  }
beta = data.matrix(c(1,-1,0.5,-0.5,1, rep(0, p-q)))
theta<-matrix(1:25,nrow=p,ncol=p)
theta[lower.tri(theta)] = t(theta)[lower.tri(theta)]
alpha = data.matrix(c(1,-1,0.5,-0.5,1, rep(0, p-q)))
x<- mvrnorm(n=n,mu=mu, Sigma=Sigma)  #X~N(0,sigma^2)
x %<>% set_colnames(str_c("x",seq(ncol(x))))
sigma.square<-var(x%*%beta)/r.square
x = scale(x, T, T)
for(i in 1:p){
  for(j in 1:p){
    if(i!=j){y <- x%*%beta+1/2*theta[i,j]*x[i,]*x[,j]+rnorm(n,0,sqrt(sigma.square))}
    else{y <- x%*%beta+rnorm(n,0,sqrt(sigma.square))}
   y = scale(y, T, T)
  }
}
soft.th=function(lambda,x){
  sign(x)*pmax(abs(x)-lambda,0)
}
rlasso=function(X, y,lambda,gamma=0.5){
  
  p=ncol(X)
  n=length(y)
  eps=1
  beta=array(0, dim=p)
  theta=array(0, dim=c(p,p))
  beta.old=array(0, dim=p)
  theta.old=array(0, dim=p)
  yhat = sum(X%*%beta)
  while(eps>0.001){
    for(j in 1:p){
      r1= y-X[,-j]%*%beta[-j]
      beta = gamma*soft.th(t(X[,j])%*%r1,lambda-alpha[j])+(1-gamma)*(t(X[,j])%*%r1)
      for(k in 1:p){
        if(k!=j){
          r2=y-yhat+(X[,j]*X[,k])*theta[j,k]
          theta= (gamma*soft.th(t(X[,j]*X[,k])%*%r2,lambda+alpha[j]+alpha[k])+(1-gamma)*(t(X[,j]*X[,k])%*%r2))/(sum(X[,j]*X[,k]))^2
        }
        else{theta = 0}
       }
      eps=max(abs(beta-beta.old),abs(theta-theta.old))
      beta.old=beta
      theta.old=theta
    }
    
    #beta.0= y-sum(X%*%beta)
    return(list(beta=beta,theta=theta))
  }
  
}
  
rlasso(X=x,y=y,lambda=1,gamma=0.5)



