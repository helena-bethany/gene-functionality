library(randomForest)

randomforestx <- function(file1){

data <- read.csv(file1, header=TRUE)
data[is.na(data)] <- 0
d <- Sys.Date()
  
ind <- sample(2,nrow(data),replace=TRUE,prob=c(0.7,0.3))
trainData <- data[ind==1,]
testData <- data[ind==2,]
data.rf <- randomForest(Functionality~.,data=trainData,ntree=1000,proximity=TRUE)
total <- round(importance(data.rf),4)
err <- data.rf$err.rate[1000,]
dataPred <- predict(data.rf, newdata=testData)
pred <- table(dataPred,testData$Functionality)
n <- 1

while (n < 100) {

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

# Calculate averages across 100 RF simulations for OOB
  
oob = mean(err$OOB)
no = mean(err$No)
yes = mean(err$Yes)
err_final = rbind(err, c(NA,oob,no,yes)) 
  
# Calculate averages across 100 RF simulations for variable importance
  
rownames(total) <- total[,1]
num <- dim(total)[1]
n <- 2
v = c()
  
while (n < num) {
i <- rowMeans(total[n,2:101])  
v <- c(v, as.numeric(i))
n <- n + 1  
}
  
total$Importance <- v
  
# Calculate prediction success rate
  
success <- pred[2,3] + pred[1,2]
total_preds <- success + pred[2,2] + pred[1,3]
success_rate <- success/total_preds
pred$Prediction <- c("",success_rate)
  
# Export files

name=substr(file1,8,1000)
  
file2=paste(d,"importance",name,sep='-')
file3=paste(d,"error",name,sep='-')
file4=paste(d,"predictions",name,sep='-')

write.csv(total,file2)  
write.csv(err_final,file3)    
write.csv(pred,file4)   

}
