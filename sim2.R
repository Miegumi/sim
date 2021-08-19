library(glmnet)
library(devtools)
# install_github('luyajun01/screening')
library(screening)
library(dplyr)
library(myscreening)
library(MASS)
library(doParallel)
registerDoParallel(cores=4)
n = 200 #data size
p = 10000 #data dim
r=0.3 #协方差rou
r.square<-0.5   ##信噪比
p1 = 5 #有效特征数
sim2 = function(n,p,mu,Sigma,r.square,r){
   mu<- rep(0,p)  #均值
   Sigma <- matrix(r, ncol=p,nrow=p)  #r是协方差
   diag(Sigma)<-1 # 协方差阵的对角线更正为1
   x<- mvrnorm(n=n,mu=mu, Sigma=Sigma)  # 产生服从N（0，Sigmas)的随机数
   beta= c(rep(5,p1), rep(0, p-p1))
   sigma.square<-var(x%*%beta)/r.square
   y <-x%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
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
result_holp = foreach(i=1:100, .combine = "rbind") %do% sim2(n,p,mu,Sigma,r.square,r)$lag_holp ##rerun 100 times
result_sis = foreach(i=1:100, .combine = "rbind") %do% sim2(n,p,mu,Sigma,r.square,r)$lag_sis
result_rrcs = foreach(i=1:100, .combine = "rbind") %do% sim2(n,p,mu,Sigma,r.square,r)$lag_rrcs
result_isis = foreach(i=1:100, .combine = "rbind") %do% sim2(n,p,mu,Sigma,r.square,r)$lag_isis
result_forward = foreach(i=1:100, .combine = "rbind") %do% sim2(n,p,mu,Sigma,r.square,r)$lag_forward

sum(result_holp)/100
sum(result_sis)/100
sum(result_rrcs)/100
sum(result_isis)/100
sum(result_forward)/100

