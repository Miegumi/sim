library(MASS)
install_github('wwrechard/screening')
library(screening)
 p=100;n=20;r=0.3
 lag<-numeric(100)
 for(j in 1:100){
    mu<- rep(0,p)  #均值
    Sigma <- matrix(NA, ncol=p,nrow=p)
  for(m in 1:p){
    for(n in 1:p){
       Sigma[m,n]=r^abs(m-n)
  }
}
  diag(Sigma)<-1  # 协方差阵的对角线更正为1
  X<- mvrnorm(n=n,mu=mu, Sigma=Sigma)  # 产生服从N（0，Sigmas)的随机数
  beta<-numeric(p)
  beta[1]=3;beta[4]=1.5;beta[7]=2
 r.square<-0.5
 sigma.square<-var(X%*%beta)/r.square
  y<-X%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
  output= screening(X, y, method = 'holp', num.select = n, ebic = TRUE)$screen
  lag[j]<-sum(c(1,4,7) %in% output)==3
  }
 sum(lag)/100
