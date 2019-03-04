randomforestx <- function(file1,file2,file3,file4,x){

library(randomForest)
data <- read.csv(file1, header=TRUE)
ind <- sample(2,nrow(data),replace=TRUE,prob=c(0.7,0.3))
trainData <- data[ind==1,]
testData <- data[ind==2,]
data.rf <- randomForest(Functionality~.,data=trainData,ntree=1000,proximity=TRUE)
total <- round(importance(data.rf),4)
err <- data.rf$err.rate[1000,]
dataPred <- predict(data.rf, newdata=testData)
pred <- table(dataPred,testData$Functionality)
n <- 1

while (n < x) {

ind <- sample(2,nrow(data),replace=TRUE,prob=c(0.7,0.3))
trainData <- data[ind==1,]
testData <- data[ind==2,]
data.rf <- randomForest(Functionality~.,data=trainData,ntree=1000,proximity=TRUE)
value <- round(importance(data.rf),4)
total <- cbind(total,value)
err2 <- data.rf$err.rate[1000,]
err <- rbind(err,err2)
dataPred <- predict(data.rf, newdata=testData)
pred2 <- table(dataPred,testData$Functionality)
pred <- pred + pred2
n <- n + 1

}

write.csv(total,file2)
write.csv(err,file3)
write.csv(pred,file4)
}