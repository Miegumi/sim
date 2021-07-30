install_github('wwrechard/screening')
library(screening)
n = 20
p = 100
times = 100
x = matrix(NA, nrow=n,ncol=p)
s=5
z = matrix(rnorm(n*p), nrow = n,ncol = p)
w = matrix(rnorm(n*p), nrow = n,ncol = p)
for(j in 1:times){
  for(i in 1:5){
    x[,1]=(z[,1]+x[,1])/sqrt(2)
    x[,i+s]=x[,i]+rnorm(n,0,0.01)
    x[,i+2*s]=x[,i]+rnorm(n,0,0.01)
  }
  for(i in 16:p){
    for(j in 1:5){
      x[,i]=(z[,i]+sum(w[,j]))/2
    }
  }
beta<-numeric(p)
for(i in 1:5){
  beta[i]<-5
}
r.square<-0.5
sigma.square<-var(x%*%beta)/r.square
y <-X%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
output= screening(X, y, method = 'holp', num.select = n, ebic = TRUE)$screen
lag[j]<-sum(c(1,2,3,4,5) %in% output)==5
}
sum(lag)/100
  
