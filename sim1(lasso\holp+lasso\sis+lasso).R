library(glmnet)
library(screening)
library(dplyr)
library(myscreening)
library(doParallel)
library(lyjtools)
registerDoParallel(cores=4)
n=200  #data size
p=10000 #data dim
p1 = 5 #有效特征数
r.square<-0.5
beta<-numeric(p)
u<-rbinom(p1,1,0.4)   #u服从伯努利分布
for(i in 1:p1){
  beta[i]<-(-1)^u[i]*(abs(rnorm(1))+4*log(n)/sqrt(n))   #生成beta
}
#sim1 = function(n,p,p1,r.square){
  x=matrix(NA,ncol=p,nrow=n)
  for (j in 1:p){
    x[,j]=rnorm(n,0,1)
  }
  sigma.square<-var(x%*%beta)/r.square
  y <- x%*%beta+rnorm(n,0,sqrt(sigma.square))   #拟合回归方程
  output_holp = myscreening::my_screening(x, y, method = 'holp', num.select = n, ebic = TRUE)$screen #这是r 方法
  output_sis = myscreening::my_screening(x, y, method = 'sis', num.select = n, ebic = TRUE)$screen
  train<-sample(1:n,round(n*0.7)) #70%为训练集,30%为测试集
  #lasso
    x_train=as.matrix(x[train,])
    x_test=as.matrix(x[-train,])
    y_train=as.numeric(y[train])
    y_test=as.numeric(y[-train])
    cv.fit<-cv.glmnet(x_train,y_train,alpha=1,family="gaussian")
    test_pre<- predict(cv.fit,newx=x_test,s=cv.fit$lambda.min)
    mse=sum((test_pre-y_test)^2)/length(y_test) 
  
  #holp+lasso
   x_holp=x[,output_holp]  #提取holp筛选出的x
   x_holp_train=as.matrix(x_holp[train,])
   x_holp_test=as.matrix(x_holp[-train,])
   cv.fit_holp<-cv.glmnet(x_holp_train,y_train,alpha=1,family="gaussian")
  test_pre_holp<- predict(cv.fit_holp,newx=x_holp_test,s=cv.fit_holp$lambda.min)
    mse_holp=sum((test_pre_holp-y_test)^2)/length(y_test)
    
  #sis+lasso
    x_sis=x[,output_sis]  #提取holp筛选出的x
    x_sis_train=as.matrix(x_sis[train,])
    x_sis_test=as.matrix(x_sis[-train,])
    cv.fit_sis<-cv.glmnet(x_sis_train,y_train,alpha=1,family="gaussian")
    test_pre_sis<- predict(cv.fit_sis,newx=x_sis_test,s=cv.fit_sis$lambda.min)
    mse_sis=sum((test_pre_sis-y_test)^2)/length(y_test)
    
    
    return(list(mse=mse,mse_holp=mse_holp,mse_sis=mse_sis))
}

result= foreach(i=1:100, .combine = "rbind") %do% sim1(n,p,p1,r.square) ##rerun 100 times
sum(result$mse)/100
sum(result$mse_holp)/100
sum(result$mse_sis)/100
