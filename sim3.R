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
r.square<-0.5  ##信噪比
sim3 = function(n,p,mu,Sigma,r.square,r){
   mu<- rep(0,p)  #均值
   Sigma <- matrix(NA, ncol=p,nrow=p)
   for(m in 1:p){
     for(n in 1:p){
       Sigma[m,n]=r^abs(m-n)   
     }
   }
   diag(Sigma)<-1  # 协方差阵的对角线更正为1
   x<- mvrnorm(n=n,mu=mu, Sigma=Sigma)  # 产生服从N（0，Sigmas)的随机数
   beta<-numeric(p)
   beta[1]=3;beta[4]=1.5;beta[7]=2
   
   sigma.square<-var(x%*%beta)/r.square
   y<-X%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
   output= screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen
   lag <-sum(c(1,4,7) %in% output)==3
   return(lag)
 }
 result = foreach(i=1:100, .combine = "rbind") %do% sim3(n,p,mu,Sigma,r.square,r) ##rerun 100 times
 sum(result)/100
