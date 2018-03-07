# Project: best design for Q-matrix

# Generate data
require(dina)
N = 400
K = 4
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
CLs = rep(1:(2^K),25)
Profiles = t(matrix(rep(t(As), 100), 4,400))
CLs.onehot <- t(matrix(rep(diag(1,16), 25), 16,400))

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
Q.1 <- Q.complete[rep(1:4,8),]
Q.2 <- Q.complete[rep(c(1:4,15),6),]
Q.3 <- Q.complete[rep(1:15,2),]


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
write(toWrite1, "LossStrategy1Skill4.txt")

toWrite2 <- sapply(c(0.01,0.05,0.1,0.2,0.3), function(i){
  myResult <- component(Q.2, CLs.onehot, As, Profiles, theta = i)
  return(colSums(myResult)/N)
})
write(toWrite2, "LossStrategy2Skill4.txt")

toWrite3 <- sapply(c(0.01,0.05,0.1,0.2,0.3), function(i){
  myResult <- component(Q.3, CLs.onehot, As, Profiles, theta = i)
  return(colSums(myResult)/N)
})
write(toWrite3, "LossStrategy3Skill4.txt")

# Draw Graph for experiment 1
t1 <- matrix(scan("LossStrategy1Skill4.txt"), 32, 5)[1:30,]
t2 <- matrix(scan("LossStrategy2Skill4.txt"), 30, 5)
t3 <- matrix(scan("LossStrategy3Skill4.txt"), 30, 5)

t1 <- cbind(1,c(1:30),t1)
t2 <- cbind(2,c(1:30),t2)
t3 <- cbind(3,c(1:30),t3)
tAll <- rbind(t1,t2,t3)
colnames(tAll) <- c("strategy", "configuration", "0.01", "0.05", "0.1", "0.2", "0.3")

library(tidyr)
d <- gather(as.data.frame(tAll), slipAndGuess, value, -strategy, -configuration)
d$strategy <- as.factor(d$strategy)

library(ggplot2)
ggplot(d, aes(configuration,value, shape=slipAndGuess, col = strategy)) + 
  geom_point()+geom_line()+labs(x = "number of questions", y = "loss")



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
theta = 0.2
J = 5
ss = gs = rep(theta,J)
start <- Sys.time()
myResult <- getAllAccuracy(J,15)
end <- Sys.time()
end - start
write(myResult, paste0("Skill4LossBestConfig_theta",theta,"_J",J,".txt"))

myResult <- matrix(scan("Skill4LossBestConfig_theta0.2_J5.txt"),6,11628)

allIndex <- getIndex(5,15)
myFactor <- factor(rep(6,ncol(myResult)), levels = c(1:6))
myIndex1 <- apply(allIndex,2,function(m) all(c(1,2,3,4) %in% m))
myFactor[myIndex1] <- 1
myIndex2 <- apply(allIndex,2,function(m) {all(c(1,2,3,4) %in% m) & any(c(5,6,7,8,9,10,11,12,13,14,15) %in% m)})
myFactor[myIndex2] <- 2
myIndex3 <- apply(allIndex,2,function(m) length(table(m))==1)
myFactor[myIndex3] <- 3
myIndex4 <- apply(allIndex,2,function(m) length(table(m))==2)
myFactor[myIndex4] <- 4
myIndex5 <- apply(allIndex,2,function(m) length(table(m))==3)
myFactor[myIndex5] <- 5

myData <- data.frame(configuration = 1:ncol(myResult), type = myFactor, loss = myResult[1,])

ggplot(myData, aes(configuration, loss, color = type)) + geom_point() +
  scale_color_manual(labels = c("I", "II","III","IV","V", "VI"), values = c("red", "blue", "black", "purple", "turquoise", "green"))
