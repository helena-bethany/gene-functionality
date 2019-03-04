#Copying parts of Paul's Script = not so fun times

#library(MASS)
#library(RColorBrewer)
#library(fields)
#library(gplots)

#k <- 11
#my.cols <- rev(brewer.pal(k, "RdYlBu"))
#nulldata <- read.csv("190131nullRF1.csv", header=TRUE)
#regM <- lm(phyloPMean ~ PhastConsMean + phyloPMax + PhastConsMax + RIsearchMIN + X.RIsearchAVE + Fu.Li.D + Fu.Li.F + SNPsave1000 + SNPs1000 + Tajima.D + aveMAF + minMAF + TTratio + transversion + transition + CM-HMM + HMM + CMx + effnseq + nseq + averageCov + incompatiable + compatiable + sigpair + covMax10 + covMin10 + AveragebpTYPE + AverageTYPE + CardiacMuscle + HematopoieticMultipotent + Myotube + BipolarNeuron + Myocyte + NeuralProgenitorCell + Hepatocyte + SmoothMuscleCell + AveragebpLINE + AverageLINE + GM12878 + K562 + hESC + HepG2 + Repeats100Query + RepeatsChr + RepeatsTotal + CM, data=nulldata)
#summary(regM)


allData <- read.csv("FinalDATA/190207nullRF1.csv", header=TRUE)
allData$Start <- NULL
allData[1:51] <- lapply(allData[1:51], as.numeric)

cormat <- round(cor(allData, method = "spearman"),2)

  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }

  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  } 

library(reshape2)

upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 9, hjust = 1))+ coord_fixed()+xlab("")+ylab("")


+annotate("segment", x=5.5, xend=5.5, y=0, yend=5.5, color="black") +  annotate("segment", x=2.5, xend=2.5, y=0, yend=2.5, color="black") + annotate("segment", x=21.5, xend=21.5, y=0, yend=21.5, color="black")
