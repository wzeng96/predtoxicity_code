## Run Functions for lasso, ridge and naive Bayesian for Tox21
# Wenyu Zeng 1/5/2023


source("runModel_functions.r")

library(caret)
library(glmnet)
library(klaR)
library(pROC)


# remove index (0,1,2,..)
data <- unique(data[,-1])


# filter data
fea <- data[,-c(1,2)]
outcome <- data[, 1]
df <- fea_to_remove(fea, 0.95)

# split data
df$outcome <- as.factor(outcome)
train <- split_data(fulldata = df, 123, prob = c(0.8,0.2))[[1]]
test <- split_data(fulldata = df, 123, prob = c(0.8,0.2))[[2]]

# Define X_train, y_train, X_test, y_test
x_train <- scale(train[,!names(train) %in% c("outcome")])
y_train <- train$outcome
x_test <- scale(test[,!names(test) %in% c("outcome")])
y_test <- test$outcome


# calculate class weights
class_weights <- weights(x_train, y_train)

# calculate lambdas
x_train_mat <- as.matrix(x_train)
y_train_mat <- as.matrix(y_train)
lambdas <- cal_lambdas(x_train_mat, y_train, weights = class_weights[[3]])
lambda_lasso <- lambdas[[1]]
lambda_ridge <- lambdas[[2]]


# Run models
weight_lasso_model <- glmnet(x_train, y_train, alpha = 1, family = quasibinomial, 
                             lambda = lambda_lasso, standardize = TRUE,
                             weights= class_weights[[3]])

weight_ridge_model <- glmnet(x_train, y_train, alpha = 0, family = quasibinomial, 
                             lambda = lambda_ridge,
                             weights= class_weights[[3]])

nb_prior <- NaiveBayes(x_train, as.factor(y_train), 
                       prior = c(class_weights[[1]], class_weights[[2]]), usekernel=TRUE)


# model predictions
options(warn=-1)
x_test <- as.matrix(x_test)
mod_pred <- model_predicts(x_test, y_test, weight_lasso_model,weight_ridge_model,nb_prior, lambdas = lambdas)

# Print lambdas
paste('Lasso lambda', lambda_lasso)
paste('Ridge Lambda', lambda_ridge)

