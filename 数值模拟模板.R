## env
library(glinternet)
library(dplyr)
library(hierNet)
library(glmnet)
library(foreach)
library(xgboost)
library(doParallel)
registerDoParallel(cores=4)

source("/home/tonylu/R/project/lib/funcs.R")

## load data
n = 3000 #data size
p = 100 #data dim
times = 100 #rerun times
## func
rerun_hier = function(n,p){
  ####################################
  ##input- n:data size p:data dimension
  ##output-model auc
  #####################################
  ##data生成
  x = matrix(rnorm(n * p), ncol = 10)
  x = scale(x, TRUE, TRUE)
  y = x[, 1] + 2 * x[, 2] + x[, 1] * x[, 2] + 3 * rnorm(n)
  y = 1 * (y > 0)
  id = seq(1, y %>% length(), 1)
  ## train_test_set split
  test_set = data.frame(id, y, x) %>%
    sample_frac(0.3)

  train_set = data.frame(id, y, x) %>%
    anti_join(test_set, by = "id")

  ## hierNet model fit
  fit_path = hierNet.logistic.path(
    x = data.matrix(train_set %>%
                      select(-y, -id)),
    y = train_set %>% pull(y)
  )

  fitcv = hierNet.cv(fit_path,
                      x = data.matrix(train_set %>%
                                        select(-y, -id)),
                      y = train_set %>% pull(y)
  )

  fit = hierNet.logistic(
    x = data.matrix(train_set %>%
                      select(-y, -id)),
    y = train_set %>% pull(y),
    lam = fitcv$lamhat.1se
  )

  ## lasso model fit
  cv_fits = cv.glmnet(
    x = data.matrix(train_set %>%
                      select(-y, -id)),
    y = train_set %>% pull(y), family = "binomial", nfolds = 5
  )

  ##predict
  yhat_hiernet <- predict(fit, data.matrix(test_set %>% select(-y, -id)))
  yhat_lasso <- predict(cv_fits, data.matrix(test_set %>% select(-y, -id)), s = cv_fits$lambda.min)

  ##compare auc
  auc_hiernet = as.numeric(aucAndKS(yhat_hiernet$prob, test_set$y)$res)
  auc_lasso = as.numeric(aucAndKS(inv_log(yhat_lasso) %>% as.vector(),test_set$y)$res)
  ##计算正确识别出交互效应次数
  inter_times = ifelse(fit$th[2,1] != 0,1,0)

  ##output
  return(list(auc_hiernet, auc_lasso, inter_times))
}

## simulation result
result = foreach(i=1:times, .combine = "rbind") %do% rerun_hier(n,p) ##rerun 100 times