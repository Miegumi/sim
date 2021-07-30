library(glmnet)
library(devtools)
#install_github('wwrechard/screening')
library(screening)
library(doParallel)
registerDoParallel(cores=4)

n=200  #data size
p=10000 #data dim
p1 = 5 #有效特征数
r.square<-0.5
sim1 = function(n,p,u,r.square){
 X=matrix(NA,ncol=p,nrow=n)
 beta<-numeric(p)
 u<-rbinom(p1,1,0.4)   #u服从伯努利分布
 for(i in 1:p1){
   beta[i]<-(-1)^u[i]*(abs(rnorm(1))+4*log(n)/sqrt(n))   #生成beta
  }
 for (i in 1:n){
   X[i,]=rnorm(p,0,1)
  }
 sigma.square<-var(X%*%beta)/r.square
 y <- X%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
 output = screening(X, y, method = 'holp', num.select = n, )$screen
 lag<-sum(seq(1,p1) %in% output)==5
 return(lag)
 }
result = foreach(i=1:100, .combine = "rbind") %do% sim1(n,p,u,r.square) ##rerun 100 times

sum(result)/100
