# Install packages
pacman::p_load('ggplot2', 'randomForest', 'tree', 'caret', 'tidyverse', 'yardstick')

# Set random seed
set.seed(1234)               

# Read in the dataset
compdat <- read.csv("data/Blastomyces_env_dataset.csv")

# Split into test and training sets
sample.ind <- sample(2, nrow(compdat), replace=T, prob = c(2/3, 1/3)) #splitting dataset
train <- compdat[sample.ind==1,] #training set
dim(train) #101 observations in training set

test<-compdat[sample.ind==2,] #test set
dim(test) # 37 observations in test
 
# Generalized linear model 
##############################################################################
model <- as.formula(paste(colnames(train)[1],"~",paste(colnames(train)[2:dim(train)[2]],collapse="+")))

glm.fit <- glm(model, family = binomial(link="logit"), data = train)
summary.aov(aov(glm.fit))

y.hat <- predict(glm.fit, newdata=test[,2:ncol(test)], type = "response")
y.hat.test <- ifelse(y.hat>0.5, 1, 0)
y.test<-test[,1]

#Mean squared error
mean((y.hat.test-y.test)^2)

# Classification tree
##############################################################################
train$BAD1_DNA_present <- factor(train$BAD1_DNA_present)
test$BAD1_DNA_present <- factor(test$BAD1_DNA_present)

tr <- tree(model,data=train)
summary(tr)

#testing error
y.hat <- predict(tr, newdata = test[,2:ncol(test)], type="class")
y.hat.test <- as.numeric(as.character(y.hat))

y.test <- test[,1]
y.test <- as.numeric(as.character(y.test))

# Calculate misclassification error
mean((y.hat.test-y.test)^2)

# Plot the calssification tree
plot(tr)
text(tr,pretty=0)

# randomForest
##############################################################################
rf <- randomForest(model, data=train, 
                   ntree=2000,
                 importance=TRUE)

plot(rf)
legend("topright", colnames(rf$err.rate), col=1:4, cex=0.8, fill=1:4)

#Plot of variable importance
varImpPlot(rf)

# Calculate training error
yhat.train <- as.numeric(as.character(rf$predicted))
y.train <- as.numeric(as.character(train[,1]))

# Mean-squared training error
mean((yhat.train-y.train)^2) 

# Tuning the model
tuneRF(x = train[,c(2:ncol(train))], y=train[,1], plot = T,
       ntreeTry = 900, improve=0.001, trace=T)


# Tuned model
rf2 <- randomForest(model, data=train,
                  ntree=2000, mtry=2,
                  importance=TRUE)


# Testing the tuned model with holdout test set
yhat.test<-as.numeric(as.character(predict(rf2,newdata=test[,2:ncol(test)])))
yhat.test<-ifelse(yhat.test > 0.5, 1, 0)
y.test <- as.numeric(as.character(test[,1]))

# Calculate testing mean-square error
mean((yhat.test-y.test)^2)  #24.3% MSE

# Confusion matrix
cm_rf <- caret::confusionMatrix(as.factor(yhat.test), as.factor(y.test))
cm_rf # 75.7% accuracy

# ROC curve
rfROC <- pROC::roc(response = yhat.test,
                    predictor = y.test)
# AUC of 0.757

plot(rfROC, type = "s", 
     col = "#529DCF", xaxs = "i", yaxs="i",
     print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.6)
