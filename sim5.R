install_github('luyajun01/screening') #安装新包
library(glmnet)
library(devtools)
library(screening)
library(doParallel)
library(dplyr)
library(myscreening)
registerDoParallel(cores=4)
n = 200 #data size
p = 10000 #data dim
p1 = 15 #有效特征数
delte.square=0.1
r.square = 0.5
sim5 = function(n,p, p1,delte.square, r.square){
  x = matrix(NA, nrow=n,ncol=p)
  z = matrix(rnorm(n*p), nrow = n,ncol = p)    ##\z_{1} ~ N(0,1)
  for(m in 0:4){
    x[,1+3*m] = z[,1]+rnorm(1,0,sqrt(delte.square))
    x[,2+3*m] = z[,2]+rnorm(1,0,sqrt(delte.square))
    x[,3+3*m] = z[,3]+rnorm(1,0,sqrt(delte.square))
  }
  for(i in (p1+1):p){
    x[,i]=rnorm(n)
  }
  x<-scale(x)
  beta = data.matrix(c(rep(3,p1), rep(0, p-p1)))
  sigma.square <- var(x%*%beta)/r.square
  y <- x%*%beta+rnorm(n,0,sqrt(sigma.square)) ##拟合方程
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

  prob_holp=length(output_holp[which(output_holp<=15)])/p1
  prob_sis=length(output_holp[which(output_sis<=15)])/p1
  prob_rrcs=length(output_holp[which(output_rrcs<=15)])/p1
  prob_isis=length(output_holp[which(output_isis<=15)])/p1
  prob_forward=length(output_holp[which(output_forward<=15)])/p1

return(list(prob_holp,prob_sis,prob_rrcs,prob_isis,prob_forward))
}
result_holp = foreach(i=1:100, .combine = "rbind") %do% sim5(n,p, p1,delte.square, r.square)$prob_holp   ##rerun 100 times
result_sis = foreach(i=1:100, .combine = "rbind") %do% sim5(n,p, p1,delte.square, r.square)$prob_sis
result_rrcs = foreach(i=1:100, .combine = "rbind") %do% sim5(n,p, p1,delte.square, r.square)$prob_rrcs
result_isis = foreach(i=1:100, .combine = "rbind") %do% sim5(n,p, p1,delte.square, r.square)$prob_isis
result_forward = foreach(i=1:100, .combine = "rbind") %do% sim5(n,p, p1,delte.square, r.square)$prob_forward

sum(result_holp)/100
sum(result_sis)/100
sum(result_rrcs)/100
sum(result_isis)/100
sum(result_forward)/100

  
 
