library(glmnet)
library(devtools)
# install_github('luyajun01/screening')
library(screening)
library(doParallel)
registerDoParallel(cores=4)
n = 200  #data size
p = 10000 #data dim
x = matrix(NA, nrow=n,ncol=p)
k = 2
p1 = 5 #有效特征数
r.square<-0.5  #信噪比
sim4 = function(n,p,k,p1,r.square){
  phi = matrix(rnorm(n*k), nrow = n,ncol = k)  ##\phi_{1} ~ N(0,1)
  f = matrix(rnorm(p*k), nrow = p,ncol = k)     ##\f_{1} ~ N(0,1)
  eta = matrix(rnorm(n*p), nrow = n,ncol = p)
  for(i in 1:p){
    for(j in 1:k){
      x[,i] = sum(phi[,j]*f[i,]) + eta[,i]
    }
  }
  
  beta = c(rep(5,5), rep(0, p-5))
  sigma.square<-var(x%*%beta)/r.square
  y <-x%*%beta+rnorm(n,0,sqrt(sigma.square)) ##拟合方程
  #output = screening::screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是r 方法
  output = myscreening::my_screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是C++ 方法
  lag <-sum(seq(1,p1) %in% output) == p1
  return(lag)
}

result = foreach(i=1:100, .combine = "rbind") %do% sim4(n,p,k,p1,r.square) ##rerun 100 times

sum(result)/100

