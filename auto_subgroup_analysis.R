# Script which divides data into training and test sets, selects a candidate 
# subgroup of interest using the training set, and tests for subgroup
# heterogeneity in the treatment effect using the test set

#setwd("~/sfg-nejm-sprint")
baseline <- read.csv("data/baseline.csv")
safety <- read.csv("data/safety.csv")
outcomes <- read.csv("data/outcomes.csv")
merged.data <- merge.data.frame(baseline, safety, by="MASKID")
merged.data <- merge.data.frame(merged.data, outcomes, by="MASKID")

outcome <- merged.data$EVENT_PRIMARY & merged.data$T_PRIMARY <= 365.25*5

#
# Auto subgroup analysis
#

# Divide dataset in half in stratified manner
library(caret)
seed = 12345; set.seed(seed)
train.ids <- createDataPartition(as.factor(outcome), p = .5, list = FALSE)
X.train <- merged.data[train.ids,]
# Map outcomes to {1,-1}
y.train <- 2*outcome[train.ids] - 1
t.train <- X.train$INTENSIVE
X.test <- merged.data[-train.ids,]

# Compute imputation values using training set
library(imputeMissings)
imputations = compute(X.train, method = "median/mode")
# Impute missing predictors in training set and in testing set
X.train.imputed = impute(X.train, imputations)
X.test.imputed = impute(X.test, imputations)

# Construct classifier g on training data to optimize
#    sum_{i=1}^n y_i  t_i g(x_i) / (sum_{j=1}^n t_j g.est(x_j) )
#  - sum_{i=1}^n y_i  (1 - t_i) g(x_i) / (sum_{j=1}^n (1-t_j) g.est(x_j) )
#  - sum_{i=1}^n y_i  t_i (1 - g(x_i)) / (sum_{j=1}^n t_j (1 - g.est(x_j) ) )
#  + sum_{i=1}^n y_i  (1-t_i) (1 - g(x_i)) / (sum_{j=1}^n (1-t_j) (1 - g.est(x_j) )
# where g.est(x_i) = indic{y_i > 0} t_i + indic{y_i < 0} (1-t_i)
g.est = (y.train > 0)*t.train + (y.train < 0)*(1-t.train) + (y.train == 0)*0.5
# Compute denominators
d11 = sum(t.train*g.est)
d01 = sum((1-t.train)*g.est)
d10 = sum(t.train*(1-g.est))
d00 = sum((1-t.train)*(1-g.est))
# Compute weighted classification parameters
z.train = y.train*t.train*(1/d11 + 1/d10) - y.train*(1-t.train)*(1/d01 + 1/d00) 

# Train xgboost classifier
library(xgboost)
variables <- c("BMI")
dtrain <- xgb.DMatrix(data = as.matrix(X.train.imputed[z.train !=0,variables]), 
                      label = (sign(z.train[z.train !=0])+1)/2, 
                      weight = abs(z.train[z.train !=0])/min(abs(z.train[z.train !=0])))
# Use CV to select number of training rounds
seed = 12345; set.seed(seed)
xgb.cv.model <- xgb.cv(data = dtrain, objective = "binary:logistic", nrounds = 10, nfold = 10)
best.nrounds = which.min(xgb.cv.model$evaluation_log$test_error_mean)
xgb.cv.model$evaluation_log[best.nrounds]
# Training final classifier using selected number of training rounds
seed = 12345; set.seed(seed)
xgb.model <- xgboost(data = dtrain, objective = "binary:logistic", nrounds = best.nrounds)
xgb.plot.tree(model = xgb.model)

# Use classifier to assign groups to test set
g.test.xgb <- predict(xgb.model, as.matrix(X.test.imputed[,variables])) >= 0.5

# Evaluate significance of interaction between group variable and treatment on test set
lm.test.model <- lm(outcome[-train.ids] ~ INTENSIVE*g.test.xgb, data=X.test)
summary(lm.test.model)
library("survival")
cox.test.model <- coxph(Surv(X.test$T_PRIMARY, outcome[-train.ids]) ~ INTENSIVE*g.test.xgb, data=X.test)
summary(cox.test.model)            

