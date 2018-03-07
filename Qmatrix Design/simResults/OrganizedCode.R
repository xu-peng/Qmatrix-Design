# Project: best design for Q-matrix

# Generate data
require(dina)
N = 200
K = 3
delta0 = rep(1,2^K)

# Defining matrix of possible attribute profiles
# As is all possible profile patterns, for 3 skills, there are 8 patterns
As = rep(0,K)
for(j in 1:K){
  temp = combn(1:K,m=j)
  tempmat = matrix(0,ncol(temp),K)
  for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
  As = rbind(As,tempmat)
}
As = as.matrix(As)

# Generate population

# Random Generation
# set.seed(1)
# PIs = rep(1/(2^K),2^K)
# CLs = c((1:(2^K)) %*% rmultinom(n=N,size=1,prob=PIs) )
# Profiles = As[CLs,]

# Deterministic Generation
CLs = rep(1:8,25)
Profiles = t(matrix(rep(t(As), 125), 3,200))
CLs.onehot <- t(matrix(rep(diag(1,8), 125), 8,200))

# Experiments 1: Three strategies
# Q.complete is consisted of all possible q-vectors, for 3 skills, there all 7 q-vectors
Q.complete = matrix(rep(diag(K),1),1*K,K,byrow=TRUE)
for(mm in 2:K){
  temp = combn(1:K,m=mm)
  tempmat = matrix(0,ncol(temp),K)
  for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
  Q.complete = rbind(Q.complete,tempmat)
}

# Three strategies
# Strategy 1: Repetition of e_{j}
# Strategy 2: Repetition of e_{j}+{1,1,1}
# Strategy 3: Repetition of the completed Q-matrix
Q.1 <- Q.complete[rep(1:3,7),]
Q.2 <- Q.complete[rep(c(1:3,7),6),]
Q.3 <- Q.complete[rep(1:7,3),]


# Function component() 
component <- function(Q.complete, CLs.onehot, As, Profiles, theta){
  result <- matrix(0, 100, nrow(Q.complete))
  for (i in 1:100){
    for (J in 1:nrow(Q.complete)){
      Q <- Q.complete[1:J,]
      if(!is.matrix(Q)){
        Q <- t(as.matrix(Q))
      }
      # Setting item parameters and generating attribute profiles
      ss = gs = rep(theta,J)
      
      # Simulate data under DINA model
      gen = DINAsim(Profiles,Q,ss,gs)
      
      
      # Calculate logLikelihood
      
      toCompare <- As %*% t(Q)
      
      correct <- t(matrix(rowSums(Q*Q),J,nrow(As)))

#     correct <- t(matrix(rep(apply(Q,1,function(x) sum(x^2)), nrow(As)), nrow(Q), nrow(As)))

      idealResponse <- apply(toCompare == correct,2,as.integer)
      
      I1 <- ((1-ss)^t(idealResponse) * (gs)^t(1-idealResponse))
      I2 <- (ss^t(idealResponse) * (1-gs)^t(1-idealResponse))
      
      logLikelihood <- t(apply(gen$Y, 1, function(m) colSums(m*log(I1) + (1-m)*log(I2))))
      toProbability <- exp(logLikelihood)
      # Since prior is uniform, we directly calulate the posterior
      posterior <- toProbability/rowSums(toProbability)
      result[i,J] <- sum((posterior - CLs.onehot)^2)
      
      # MLE <- apply(logLikelihood,1,which.max)
      # result[i,J] <- sum(MLE == CLs)
      
    }
  }
  return(result)
}

# Experiment1 Implementation
set.seed(100)
toWrite1 <- sapply(c(0.01,0.05,0.1,0.2,0.3), function(i){
              myResult <- component(Q.1, CLs.onehot, As, Profiles, theta = i)
              return(colSums(myResult)/N)
              })
write(toWrite1, "LossStrategy1.txt")

toWrite2 <- sapply(c(0.01,0.05,0.1,0.2,0.3), function(i){
  myResult <- component(Q.2, CLs.onehot, As, Profiles, theta = i)
  return(colSums(myResult)/N)
})
write(toWrite2, "LossStrategy2.txt")

toWrite3 <- sapply(c(0.01,0.05,0.1,0.2,0.3), function(i){
  myResult <- component(Q.3, CLs.onehot, As, Profiles, theta = i)
  return(colSums(myResult)/N)
})
write(toWrite3, "LossStrategy3.txt")

# Draw Graph for experiment 1
t1 <- matrix(scan("LossStrategy1.txt"), 21, 5)
t2 <- matrix(scan("LossStrategy2.txt"), 24, 5)[1:21,]
t3 <- matrix(scan("LossStrategy3.txt"), 21, 5)

t1 <- cbind(1,c(1:21),t1)
t2 <- cbind(2,c(1:21),t2)
t3 <- cbind(3,c(1:21),t3)
tAll <- rbind(t1,t2,t3)
colnames(tAll) <- c("strategy", "configuration", "0.01", "0.05", "0.1", "0.2", "0.3")

library(tidyr)
d <- gather(as.data.frame(tAll), slipAndGuess, value, -strategy, -configuration)

library(tikzDevice)
options(tz="CA")
require(ggplot2)
d$strategy <- as.factor(d$strategy)

g <- ggplot(d, aes(configuration,value, shape=slipAndGuess, col = strategy)) + 
  geom_point()+geom_line()+labs(x = "number of questions", y = "loss")

tikz(file = "LossStrategyComparison.tex", width = 4, height = 3)
print(g)
dev.off()

# n number of questions
# k number of patterns
# function getIndex() returns a matrix containing the indices for q-vectors for each pattern combination
getIndex <- function(n,k){
  # Generate all possible choice of pattern combination
  allPatterns<- combn(n+k-1,k-1)
  result <- apply(allPatterns, 2, function(m){
    temp <- c(0,m,n+k)
    pattCounts <- diff(temp)-1  #pattCounts stores the frequency for each pattern corresponding to a chosen pattern combination
    index <- Reduce(c,mapply(rep, 1:k, pattCounts))
  })
  return(result)
}

getAllAccuracy<- function(n,k){
  allIndex <- getIndex(n,k)
  # allAccurcy <- apply(allIndex, 2, getAccuracy)
  allAccuracy <- foreach(a = split(allIndex, col(allIndex)), 
                         .combine = cbind, 
                         .export=c('getAccuracy',
                                   'Q.complete',
                                   'As',
                                   'ss',
                                   'gs',
                                   'getSingleAccuracy',
                                   'DINAsim',
                                   'Profiles',
                                   'CLs.onehot')) %dopar% {getAccuracy(a)}
  allAccuracy
}

getSingleAccuracy <- function(i,Q){
  # Simulate data under DINA model
  gen = DINAsim(Profiles,Q,ss,gs)
  
  logLikelihood <- t(apply(gen$Y, 1, function(m) colSums(m*log(I1) + (1-m)*log(I2))))
  # Maximum Likelihood Estimation
  # MLE <- apply(logLikelihood,1,which.max)
  # return(sum(MLE==CLs))
  
  toProbability <- exp(logLikelihood)
  # Since prior is uniform, we directly calulate the posterior
  posterior <- toProbability/rowSums(toProbability)
  return(sum((posterior - CLs.onehot)^2))
}

getAccuracy <- function(myIndex){
  Q <- Q.complete[myIndex,]
  
  # Get ideal reponse 
  toCompare <- As %*% t(Q)
  correct <- t(matrix(rep(apply(Q,1,function(x) sum(x^2)), nrow(As)), nrow(Q), nrow(As)))
  idealResponse <- apply(toCompare == correct,2,as.integer)
  
  # Calculate logLikelihood
  I1 <<- ((1-ss)^t(idealResponse) * (gs)^t(1-idealResponse))
  I2 <<- (ss^t(idealResponse) * (1-gs)^t(1-idealResponse))
  
  avg <- mean(sapply(1:100, getSingleAccuracy, Q=Q))
  return(c(avg, myIndex))
}

# Parallel computing
library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

set.seed(1)
theta = 0.01
J = 9
ss = gs = rep(theta,J)
start <- Sys.time()
myResult <- getAllAccuracy(J,7)
end <- Sys.time()
end - start
write(myResult, paste0("LossBestConfig_theta",theta,"_J",J,".txt"))



# Factoring
# We divide all Q-matrix in 4 categories
# 1 - Complete Q-matrix using only boolean unit q-vectors
# 2 - COmplete Q-matrix using at least one other q-vector
# 3 - Incomplete Q-matrix containing only 1 type q-vector
# 4 - Incomplete Q-matrix containing only 2 type q-vectors
myResult <- matrix(scan("LossBestConfig_theta0.3_J8.txt"),9,3003)
allIndex <- getIndex(8,7)
myFactor <- factor(rep(5,ncol(myResult)), levels = c(1:5))
myIndex1 <- apply(allIndex,2,function(m) all(c(1,2,3) %in% m))
myFactor[myIndex1] <- 1
myIndex2 <- apply(allIndex,2,function(m) {all(c(1,2,3) %in% m) & any(c(4,5,6,7) %in% m)})
myFactor[myIndex2] <- 2
myIndex3 <- apply(allIndex,2,function(m) length(table(m))==1)
myFactor[myIndex3] <- 3
myIndex4 <- apply(allIndex,2,function(m) length(table(m))==2)
myFactor[myIndex4] <- 4

myData <- data.frame(configuration = 1:ncol(myResult), type = myFactor, loss = myResult[1,])

ggplot(myData, aes(configuration, loss, color = type)) + geom_point() +
  scale_color_manual(labels = c("I", "II","III","IV","V"), values = c("red", "blue", "black", "purple", "turquoise"))

require(data.table)
require(reshape2)
preprocess <- function(filename, J){
  myResult <- scan(filename)
  myResult <- matrix(myResult, 24+J, length(myResult)/(24+J))  # 29 = 8(categories)*3+5(questions)
  myFinal <- data.table(configuration = c(1:ncol(myResult)), 
                        precision1 = myResult[1,]/25, 
                        recall1 = myResult[1,]/(myResult[1,]+myResult[2,]),
                        precision2 = myResult[4,]/25,
                        recall2 = myResult[4,]/(myResult[4,]+myResult[5,]),
                        precision3 = myResult[7,]/25, 
                        recall3 = myResult[7,]/(myResult[7,]+myResult[8,]),
                        precision4 = myResult[10,]/25, 
                        recall4 = myResult[10,]/(myResult[10,]+myResult[11,]),
                        precision5 = myResult[13,]/25, 
                        recall5 = myResult[13,]/(myResult[13,]+myResult[14,]),
                        precision6 = myResult[16,]/25, 
                        recall6 = myResult[16,]/(myResult[16,]+myResult[17,]),
                        precision7 = myResult[19,]/25, 
                        recall7 = myResult[19,]/(myResult[19,]+myResult[20,]),
                        precision8 = myResult[1,]/25, 
                        recall8 = myResult[22,]/(myResult[22,]+myResult[23,]))
  
  myFinal[, fscore1 := 2*precision1*recall1/(precision1+recall1)]
  myFinal[, fscore2 := 2*precision2*recall2/(precision2+recall2)]
  myFinal[, fscore3 := 2*precision3*recall3/(precision3+recall3)]
  myFinal[, fscore4 := 2*precision4*recall4/(precision4+recall4)]
  myFinal[, fscore5 := 2*precision5*recall5/(precision5+recall5)]
  myFinal[, fscore6 := 2*precision6*recall6/(precision6+recall6)]
  myFinal[, fscore7 := 2*precision7*recall7/(precision7+recall7)]
  myFinal[, fscore8 := 2*precision8*recall8/(precision8+recall8)]
  myFinal[, totalPrecision := (precision1 + precision2 + precision3 + precision4 + precision5 + precision6 + precision7 +precision8)/8]
  myFinal[, totalRecall := (recall1 + recall2 + recall3 + recall4 + recall5 + recall6 + recall7 +recall8)/8]
  myFinal[, totalFscore := (fscore1 + fscore2 + fscore3 + fscore4 + fscore5 + fscore6 + fscore7 +fscore8)/8]
  
  # ggplot(myFinal)+geom_point(data = myFinal[,.(configuration, fscore1)])+geom_point(data = myFinal[,.(configuration, fscore2)], colour='red')

  d <- melt(myFinal, id.vars="configuration")
  d
}

J <- 4
d <- preprocess("simPerPatternResult_theta0.2_J4.txt", J)


allIndex <- getIndex(J,7)
myIndex1 <- which(apply(allIndex,2,function(m) all(c(1,2,3) %in% m)))
myIndex2 <- which(apply(allIndex,2,function(m) {all(c(1,2,3) %in% m) & any(c(4,5,6,7) %in% m)}))

tikz(file = "theta02J4Total.tex", width = 4, height = 3)
#Simple plot of the dummy data using LaTeX elements
g <- ggplot(d[d[,variable %like% "total"]], aes(configuration,value, shape = variable)) + 
  geom_point() + scale_shape_manual(values=1:9) + 
  geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "total")]],col = "red") +
  geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "total")]],col = "blue")
#This line is only necessary if you want to preview the plot right after compiling
print(g)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()

tikz(file = "theta02J4Precision.tex", width = 4, height = 3)
#Simple plot of the dummy data using LaTeX elements
g <- ggplot(d[d[,variable %like% "recision"]], aes(configuration,value, shape = variable)) + 
  geom_point() + scale_shape_manual(values=1:9) + 
  geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "recision")]],col = "red") +
  geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "recision")]],col = "blue")
#This line is only necessary if you want to preview the plot right after compiling
print(g)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()

tikz(file = "theta02J4Fscore.tex", width = 4, height = 3)
#Simple plot of the dummy data using LaTeX elements
g <- ggplot(d[d[,variable %like% "score"]], aes(configuration,value, shape = variable)) + 
  geom_point() + scale_shape_manual(values=1:9) + 
  geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "score")]],col = "red") +
  geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "score")]],col = "blue")
#This line is only necessary if you want to preview the plot right after compiling
print(g)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()

tikz(file = "theta02J4Recall.tex", width = 4, height = 3)
#Simple plot of the dummy data using LaTeX elements
g <- ggplot(d[d[,variable %like% "call"]], aes(configuration,value, shape = variable)) + 
  geom_point() + scale_shape_manual(values=1:9) + 
  geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "call")]],col = "red") +
  geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "call")]],col = "blue")
#This line is only necessary if you want to preview the plot right after compiling
print(g)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()