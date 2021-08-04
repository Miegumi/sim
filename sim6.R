library(glmnet)
  library(devtools)
  install_github('luyajun01/screening')
  library(screening)
  library(doParallel)
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
      x[,i]=(z[,i]+x[,i])/sqrt(2)
      x[,i+s]=x[,i]+rnorm(n,0,sqrt(0.01))
      x[,i+2*s]=x[,i]+rnorm(n,0,sqrt(0.01))
    }
    for(i in (p1+2*s+1):p){
      for(j in 1:5){
        x[,i]=(z[,i]+sum(w[,j]))/2
      }
    }
    beta = c(rep(5,p1), rep(0, p-p1))
    sigma.square<-var(x%*%beta)/r.square
    y <-x%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
    #output = screening::screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是r 方法
    output = myscreening::my_screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是C++ 方法
    lag<-sum(seq(1,p1) %in% output)==p1
    return(lag)
  }
  result = foreach(i=1:100, .combine = "rbind") %do% sim6(n,p,p1,s,r.square) ##rerun 100 times
  
  sum(result)/100

