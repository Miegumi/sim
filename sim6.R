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

sim6 = function(n,p,p1,s,r.square){
  x = matrix(NA, nrow=n,ncol=p)
  z = matrix(rnorm(n*p), nrow = n,ncol = p)
  w = matrix(rnorm(n*p), nrow = n,ncol = p)
  for(i in 1:p1){
    x[,i]=(z[,i]+w[,i])/sqrt(2)
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

  prob_holp=length(output_holp[which(output_holp<=5)])/p1
  prob_sis=length(output_holp[which(output_sis<=5)])/p1
  prob_rrcs=length(output_holp[which(output_rrcs<=5)])/p1
  prob_isis=length(output_holp[which(output_isis<=5)])/p1
  prob_forward=length(output_holp[which(output_forward<=5)])/p1

return(list(prob_holp,prob_sis,prob_rrcs,prob_isis,prob_forward))
}
result_holp = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$prob_holp   ##rerun 100 times
result_sis = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$prob_sis
result_rrcs = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$prob_rrcs
result_isis = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$prob_isis
result_forward = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square)$prob_forward

sum(result_holp)/100
sum(result_sis)/100
sum(result_rrcs)/100
sum(result_isis)/100
sum(result_forward)/100

  
