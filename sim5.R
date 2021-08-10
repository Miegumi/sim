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
  #output = screening::screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是r 方法
  output = myscreening::my_screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是C++ 方法
  lag <-sum(seq(1,p1) %in% output) == p1
  return(lag)
}

result = foreach(i=1:100, .combine = "rbind") %do% sim5(n,p, p1,delte.square, r.square) ##rerun 100 times

sum(result)/100
