library(glmnet)
library(devtools)
# install_github('luyajun01/screening')
library(screening)
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
   y <-X%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
   #output = screening::screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是r 方法
   output = myscreening::my_screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是C++ 方法
   output= screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen
   lag <-sum(seq(1,p1) %in% output) == p1
   return(lag)
 }
result = foreach(i=1:100, .combine = "rbind") %do% sim2(n,p,mu,Sigma,r.square,r) ##rerun 100 times

sum(result)/100
