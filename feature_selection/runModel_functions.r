## Define Functions for lasso, ridge and naive Bayesian for Tox21
# Wenyu Zeng 12/2/2022

library(caret)

# CHECK
# 1) X_train and X_test have same column names
# 2) X_train and X_test are both numeric
# 3) y_train and y_test are factorize
# 4) y_train and y_test are the same level (i.e. 0s & 1s)

## remove features to clean data
fea_to_remove <- function(trainingset_no_outcome, cutoff) {
  # No Zero Variance
  nzv <- nearZeroVar(trainingset_no_outcome, uniqueCut=1)
  print('Near Zero Variance Features:')
  print(colnames(trainingset_no_outcome[,c(nzv)]))
  
  # Linear Combination
  lc <- findLinearCombos(trainingset_no_outcome)
  # print("Suggest to remove these features due to Linear Combinations: ")
  # print(lc$remove)
  # print(colnames(trainingset_no_outcome[,c(lc$remove)]))

  # High collinearity
  hc <- sort(findCorrelation(cor(trainingset_no_outcome), cutoff=cutoff))
  # print("Features that have >95% collinearity")
  # print(hc)
  # print(colnames(trainingset_no_outcome[,c(hc)]))
  
  common_fea <- intersect(intersect(nzv,lc$remove),hc)
  # print("Common elements between the three")
  # print(common_fea)
  # print(colnames(trainingset_no_outcome[,c(common_fea)]))
  
  print('Number of features removed')
  print(length(nzv)+length(lc$remove)+length(hc))
  unimp_fea <- unique(c(colnames(trainingset_no_outcome[,c(nzv)]),
                 colnames(trainingset_no_outcome[,c(hc)]),
                 colnames(trainingset_no_outcome[,c(common_fea)])))
  new_data <- trainingset_no_outcome[,-c(nzv,lc$remove,hc)]
  return(new_data)
}
### function for split data into training and testing
split_data <- function(fulldata = new_data, seednum, prob=c(0.8,0.2)){
  set.seed(seednum)
  sample <- sample(c(TRUE, FALSE), nrow(fulldata), replace=TRUE, prob=prob)
  train <- fulldata[sample,]
  test <- fulldata[!sample,]
  #train_scale <- scale(train)
  #test_scale <- scale(test)
  return(list(train, test))
}
### Function for class weights/priors
weights <- function(x_train, y_train){
  # Calculate weights for each class
  class_0 <- nrow(x_train)/(2*table(y_train)[[1]])
  class_1 <- nrow(x_train)/(2*table(y_train)[[2]])
  # create vector of weights
  weights <- ifelse(y_train == '0', class_0, class_1)
  return(list(class_0, class_1, weights))
}

### Function for lambda
cal_lambdas <- function(x_train, y_train, weights = weights[[3]]){
  lambdas <- 10^seq(2, -3, by = -.1)
  weight_model_lasso <- cv.glmnet(x_train, y_train, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5, 
                                  family = quasibinomial, weights = weights)

  # Best 
  lambda_lasso <- weight_model_lasso$lambda.min 

  
  weight_model_ridge <- cv.glmnet(x_train, y_train, alpha = 0, lambda = lambdas, standardize = TRUE, nfolds = 5, 
                                  family = quasibinomial, weights = weights)
  # Best
  lambda_ridge <- weight_model_ridge$lambda.min
  
  return(list(lambda_lasso, lambda_ridge))
}


## function for predict models on testing set
model_predicts <- function(x_test, y_test, weight_lasso_model,weight_ridge_model, nb_prior, lambdas){
  #PREDICTION
  
  # Lasso
  predictions_test <- predict(weight_lasso_model, s = lambdas[[1]], newx = x_test, type = "response")
  predicted.classes_weight <- ifelse(predictions_test > 0.5, "1", "0")
  cm <- confusionMatrix(as.factor(as.data.frame(predicted.classes_weight)$s1), y_test)
  lasso_ba <- cm$byClass[['Balanced Accuracy']]
  lasso_acc <- cm$overall[['Accuracy']]
  lasso_precision <- cm$byClass[['Precision']]
  lasso_recall <- cm$byClass[['Recall']]
  lasso_auc <- roc(y_test, predictions_test)$auc
  
  paste('Balanced Accuracy: ',lasso_ba)
  print('Lasso')
  
  print(lasso_auc)
  
  # Ridge
  predictions_test <- predict(weight_ridge_model, s = lambdas[[2]], newx = x_test, type = "response")
  predicted.classes_weight <- ifelse(predictions_test > 0.5, "1", "0")
  cm <- confusionMatrix(as.factor(as.data.frame(predicted.classes_weight)$s1), y_test)
  ridge_ba <- cm$byClass[['Balanced Accuracy']]
  ridge_acc <- cm$overall[['Accuracy']]
  ridge_precision <- cm$byClass[['Precision']]
  ridge_recall <- cm$byClass[['Recall']]
  ridge_auc <- roc(y_test, predictions_test)$auc
  
  print('Ridge')
  paste('Balanced Accuracy: ',ridge_ba)
  print(ridge_auc)
  
  # Naive Bayes with Priors
  nb_prior_pred <- predict(nb_prior, x_test)
  cm <- confusionMatrix(data = nb_prior_pred$class, reference = as.factor(y_test))
  nb_prior_ba <- cm$byClass[['Balanced Accuracy']]
  nb_prior_acc <- cm$overall[['Accuracy']]
  nb_prior_precision <- cm$byClass[['Precision']]
  nb_prior_recall <- cm$byClass[['Recall']]
  nb_prior_auc <- roc(y_test, nb_prior_pred$posterior[,2])$auc
  
  print('Naive Bayes')
  paste('Balanced Accuracy: ',nb_prior_ba)
  print(nb_prior_auc)
  
  ba <- c(lasso_ba,ridge_ba,nb_prior_ba)
  auc <- c(lasso_auc, ridge_auc, nb_prior_auc)
  acc <- c(lasso_acc, ridge_acc, nb_prior_acc)
  precision <- c(lasso_precision, ridge_precision, nb_prior_precision)
  recall <- c(lasso_recall, ridge_recall, nb_prior_recall)
  return(list(ba, auc, acc, precision, recall))
}
