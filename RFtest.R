
rftest <- function(file1,file2,file3,file4,file5,x){

library(randomForest)
data <- read.csv(file1, header=TRUE)
data2 <- read.csv(file2, header=TRUE)

ind <- sample(2,nrow(data),replace=TRUE,prob=c(0.7,0.3))
train <- data[ind==1,]

ind2 <- sample(2, nrow(data2), replace=TRUE, prob=c(0.7,0.3))
testData <- data2[ind2==2,]
train2 <- data2[ind2==1,]

trainData <- rbind(train,train2)

#trainData <- data
#testData <- data2

data.rf <- randomForest(Functionality~.,data=trainData,ntree=1000,proximity=TRUE)
dataPred <- predict(data.rf, newdata=testData)
pred <- table(dataPred,testData$Functionality)
total <- round(importance(data.rf),4)
err <- data.rf$err.rate[1000,]
n <- 1

while (n < x) {

ind <- sample(2,nrow(data),replace=TRUE,prob=c(0.7,0.3))
train <- data[ind==1,]

ind2 <- sample(2, nrow(data2), replace=TRUE, prob=c(0.7,0.3))
testData <- data2[ind2==2,]
train2 <- data2[ind2==1,]

trainData <- rbind(train,train2)

data.rf <- randomForest(Functionality~.,data=trainData,ntree=1000,proximity=TRUE)
dataPred <- predict(data.rf, newdata=testData)
pred2 <- table(dataPred,testData$Functionality)
pred <- pred + pred2
value <- round(importance(data.rf),4)
total <- cbind(total,value)
err2 <- data.rf$err.rate[1000,]
err <- rbind(err,err2)
n <- n + 1

}

write.csv(total,file3)
write.csv(err,file4)
write.csv(pred,file5)

}



