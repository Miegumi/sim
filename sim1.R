library(glmnet)
library(devtools)
install_github('wwrechard/screening')
library(screening)
p=100;n=20;
output.holp=list()
set.seed(1)
lag<-numeric(100)
for(j in 1:100){
  X=matrix(NA,ncol=p,nrow=n)
  beta<-numeric(p)
  u<-rbinom(5,1,0.4)   #u服从伯努利分布
  for(i in 1:5){
    beta[i]<-(-1)^u[i]*(abs(rnorm(1))+4*log(n)/sqrt(n))   #生成beta
  }
  
  for (i in 1:n){
    X[i,]=rnorm(p,0,1)
  }
  r.square<-0.5
  sigma.square<-var(X%*%beta)/r.square
  y <- X%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
  output.holp[[j]]= screening(X, y, method = 'rrcs', num.select = n,ebic = TRUE )
  output = screening(X, y, method = 'holp', num.select = n, )$screen
  lag[j]<-sum(c(1,2,3,4,5) %in% output)==5
}
sum(lag)/100
