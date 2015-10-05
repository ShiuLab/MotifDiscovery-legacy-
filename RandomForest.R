#!/usr/bin/env Rscript 
# Uses a presence/absence of motif(s) data frame to predict Class using RandomForest
# and determine importance.

# USAGE
# RandomForest.R [dataframe]

library("randomForest")

cArgs <- commandArgs(trailingOnly=TRUE)
file <- cArgs[1]
df <- read.table(file, header=TRUE, sep = "\t")
df$Class <- as.factor(df$Class)

"Calculating Importance..."
train_4importance <- randomForest(df[,c(3:length(df))], df$Class,ntree=500, importance = TRUE, do.trace=100)
imp <- importance(train_4importance)
savename <- sub('.txt','.imp.txt',file)
write.table(imp, savename, sep="\t")

"Calculating F measure..."
reptimes = 10
train <- replicate(reptimes, rfcv(df[,c(3:length(df))], df$Class, cv.fold=10, scale='log',step = .5, recursive=FALSE), simplify=FALSE)
predicted.cv <- sapply(train, "[[", "predicted")
ClassCol <- which(names(df)=="Class")

f<-c()
for(i in 1:reptimes){
  data <- unlist(predicted.cv[,i][1])
  TP <- sum(sapply(1:length(data), function(x) if(data[x] == "1" && df[x,ClassCol] == "1") 1 else 0))
  TN <- sum(sapply(1:length(data), function(x) if(data[x] == "0" && df[x,ClassCol] == "0") 1 else 0))
  FP <- sum(sapply(1:length(data), function(x) if(data[x] == "1" && df[x,ClassCol] == "0") 1 else 0))
  FN <- sum(sapply(1:length(data), function(x) if(data[x] == "0" && df[x,ClassCol] == "1") 1 else 0))
  p <- TP/(TP+FP)
  r <- TP/(TP+FN)
  f <- c(f,2*p*r/(p+r))
}

resultsname <- sub('.txt','.Results.txt',file)
sink(resultsname)

"Mean: "
mean(f)
"St.dev: "
sd(f)
"St.error:"
sd(f)/sqrt(reptimes)
a <- 0.975
error <- qnorm(a)*sd(f)/sqrt(reptimes)
"95% CI_lower"
mean(f)-error
"95% CI_upper"
mean(f)+error

sink()
"Done :)"

