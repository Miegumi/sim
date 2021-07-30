library(MASS)
install_github('wwrechard/screening')
library(screening)
p=1000;n=20;r=0.3
 set.seed(1)
 lag<-numeric(100)
 for(j in 1:100){
   
   mu<- rep(0,p)  #均值
   Sigma <- matrix(r, ncol=p,nrow=p) #协方差r=0.3
   diag(Sigma)<-1 # 协方差阵的对角线更正为1
   X<- mvrnorm(n=n,mu=mu, Sigma=Sigma)  # 产生服从N（0，Sigmas)的随机数
   
   beta<-numeric(p)
   for(i in 1:5){
     beta[i]<-5   #生成beta
   }
   r.square<-0.5
   sigma.square<-var(X%*%beta)/r.square
   y <-X%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
   output= screening(X, y, method = 'holp', num.select = n, ebic = TRUE)$screen
   lag[j]<-sum(c(1,2,3,4,5) %in% output)==5
 }
 sum(lag)/100
