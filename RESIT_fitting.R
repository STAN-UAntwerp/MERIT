####
# Linear Regression with factors
####
train_linear <- function(X,y,pars = list(),includeY)
{
  data <- data.frame(y,X)
  if (includeY > 0){
    for (col in ncol(data):(ncol(data) - includeY + 1)){
      data[,col] <- as.factor(data[,col])
    }
    
  }
  mod <- lm(y ~ ., data = data)
  result <- list()
  result$Yfit = as.matrix(mod$fitted.values)
  result$residuals = as.matrix(mod$residuals)
  result$model = mod
  return(result)
}

####
# random forest regression 
####
train_randomforest <- function(X,y,pars = list(),includeY = T){
  data = data.frame(y,X)
  mod <- randomForest::randomForest(y ~ ., data = data, mtry = ncol(X), ntree = 100)
  result <- list()
  result$Yfit <- as.matrix(mod$predicted)
  result$residuals <- as.matrix(y - mod$predicted)
  result$model <- mod
  return(result)
}

####
# trian model
####
train_model <- function(f,X,y,pars = list(),includeY)
{
  result <- f(X,y,pars,includeY)
}