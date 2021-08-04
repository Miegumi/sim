library(glmnet)
library(devtools)
# install_github('luyajun01/screening')
library(screening)
library(doParallel)
registerDoParallel(cores=4)
n=200  #data size
p=10000 #data dim
p1 = 5 #有效特征数
r.square<-0.5
sim1 = function(n,p,u,r.square){
 x=matrix(NA,ncol=p,nrow=n)
 beta<-numeric(p)
 u<-rbinom(p1,1,0.4)   #u服从伯努利分布
 for(i in 1:p1){
   beta[i]<-(-1)^u[i]*(abs(rnorm(1))+4*log(n)/sqrt(n))   #生成beta
  }
 for (i in 1:n){
   x[i,]=rnorm(p,0,1)
  }
 sigma.square<-var(x%*%beta)/r.square
 y <- x%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
 #output = screening::screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是r 方法
 output = myscreening::my_screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是C++ 方法
 output = screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen
 lag<-sum(seq(1,p1) %in% output)==p1
 return(lag)
 }
result = foreach(i=1:100, .combine = "rbind") %do% sim1(n,p,u,r.square) ##rerun 100 times

sum(result)/100
