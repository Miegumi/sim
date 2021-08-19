library(glmnet)
library(devtools)
library(dplyr)
library(myscreening)
## install_github('luyajun01/screening')
library(screening)
library(doParallel)
library(dplyr)
library(myscreening)
registerDoParallel(cores=4)
n = 200
p = 10000
s=5
p1 = 5 #有效特征数
r.square=0.5

##function debug 时,你可以将function 注释掉,一步一步跑看下每步结果是否正常

#sim6 = function(n,p,p1,s,r.square){
  x = matrix(NA, nrow=n,ncol=p)
  z = matrix(rnorm(n*p), nrow = n,ncol = p)
  w = matrix(rnorm(n*p), nrow = n,ncol = p)
  for(i in 1:p1){
    x[,i]=(z[,i]+x[,i])/sqrt(2)
    x[,i+s]=x[,i]+rnorm(n,0,sqrt(0.01))
    x[,i+2*s]=x[,i]+rnorm(n,0,sqrt(0.01))
  }
  for(i in (p1+2*s+1):p){
    for(j in 1:5){
      x[,i]=(z[,i]+sum(w[,j]))/2
    }
  }
  x<-scale(x)
  beta = data.matrix(c(rep(5,p1), rep(0, p-p1)))
  sigma.square<-var(x%*%beta)/r.square ##这个方差全为NA,有问题
  y <-x%*%beta + rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
  ##以下为r方法
  #output_holp = screening::screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是r 方法
  #output_sis = screening::screening(x, y, method = 'sis', num.select = n, ebic = TRUE)$screen  #这是r 方法
  #output_rrcs = screening::screening(x, y, method = 'rrcs', num.select = n, ebic = TRUE)$screen  #这是r 方法
  #output_isis = SIS::SIS(x,y,family='gaussian', tune='bic',iter=TRUE)$ix    #这是r 方法
  #output_forward = screening::screening(x, y, method = 'forward', num.select = n, ebic = TRUE)$screen  #这是r 方法
  ##以下为c++方法
  output_holp = myscreening::my_screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是C++ 方法
  output_sis = myscreening::my_screening(x, y, method = 'sis', num.select = n, ebic = TRUE)$screen #这是C++ 方法
  output_rrcs = myscreening::my_screening(x, y, method = 'rrcs', num.select = n, ebic = TRUE)$screen #这是C++ 方法
  output_isis = SIS::SIS(x,y,family='gaussian', tune='bic',iter=TRUE)$ix #这是C++ 方法
  output_forward = myscreening::my_screening(x, y, method = 'forward', num.select = n, ebic = TRUE)$screen #这是C++ 方法

  lag_holp<-sum(seq(1,p1) %in% output_holp)==p1
  lag_sis<-sum(seq(1,p1) %in% output_sis)==p1
  lag_rrcs<-sum(seq(1,p1) %in% output_rrcs)==p1
  lag_isis<-sum(seq(1,p1) %in% output_isis)==p1
  lag_forward<-sum(seq(1,p1) %in% output_forward)==p1

return(list(lag_holp,lag_sis,lag_rrcs,lag_isis,lag_forward))
}
result_holp = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$lag_holp ##rerun 100 times
result_sis = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$lag_sis
result_rrcs = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$lag_rrcs
result_isis = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$lag_isis
result_forward = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$lag_forward

sum(result_holp)/100
sum(result_sis)/100
sum(result_rrcs)/100
sum(result_isis)/100
sum(result_forward)/100
