# Project: best design for Q-matrix

# Generate data
require(dina)
N = 200
K = 3
delta0 = rep(1,2^K)

# Creating Q matrix
# Q = matrix(rep(diag(K),2),2*K,K,byrow=TRUE)
# for(mm in 2:K){
#   temp = combn(1:K,m=mm)
#   tempmat = matrix(0,ncol(temp),K)
#   for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
#   Q = rbind(Q,tempmat)
# }
# Q = Q[1:J,]

# Q.complete = diag(K);
Q.complete = matrix(rep(diag(K),1),1*K,K,byrow=TRUE)
for(mm in 2:K){
  temp = combn(1:K,m=mm)
  tempmat = matrix(0,ncol(temp),K)
  for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
  Q.complete = rbind(Q.complete,tempmat)
}

# Q <- Q.complete[c(2,3,5),]

# # Defining matrix of possible attribute profiles
As = rep(0,K)
for(j in 1:K){
  temp = combn(1:K,m=j)
  tempmat = matrix(0,ncol(temp),K)
  for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
  As = rbind(As,tempmat)
}
As = as.matrix(As)

Q.1 <- Q.complete[rep(1:3,7),]
Q.2 <- Q.complete[rep(c(1:3,7),6),]
Q.3 <- Q.complete[rep(1:7,3),]

# Generate population
set.seed(1)
PIs = rep(1/(2^K),2^K)
CLs = c((1:(2^K)) %*% rmultinom(n=N,size=1,prob=PIs) )
# Profiles = As[CLs,]
CLs = rep(1:8,25)
Profiles = t(matrix(rep(t(As), 125), 3,200))
theta = 0.1
J = 4
ss = gs = rep(theta,J)

component <- function(Q.complete, CLs, As, Profiles, theta){
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
      
      ptm <- proc.time()
      correct2 <- t(matrix(rowSums(Q*Q),J,nrow(As)))
      proc.time() - ptm
      
      ptm <- proc.time()
      correct <- t(matrix(rep(apply(Q,1,function(x) sum(x^2)), nrow(As)), nrow(Q), nrow(As)))
      proc.time() - ptm
      
      idealResponse <- apply(toCompare == correct,2,as.integer)
      
      I1 <- ((1-ss)^t(idealResponse) * (gs)^t(1-idealResponse))
      I2 <- (ss^t(idealResponse) * (1-gs)^t(1-idealResponse))
      
      logLikelihood <- t(apply(gen$Y, 1, function(m) colSums(m*log(I1) + (1-m)*log(I2))))
      
      MLE <- apply(logLikelihood,1,which.max)
      result[i,J] <- sum(MLE == CLs)
    }
  }
  return(result)
}

set.seed(100)
myResult1 <- component(Q.3, CLs, As, Profiles, theta = 0.01)
colSums(myResult1)/N
set.seed(100)
myResult2 <- component(Q.3, CLs, As, Profiles, theta = 0.05)
colSums(myResult2)/N
set.seed(100)
myResult3 <- component(Q.3, CLs, As, Profiles, theta = 0.1)
colSums(myResult3)/N
set.seed(100)
myResult4 <- component(Q.3, CLs, As, Profiles, theta = 0.2)
colSums(myResult4)/N
set.seed(100)
myResult5 <- component(Q.3, CLs, As, Profiles, theta = 0.3)
colSums(myResult5)/N

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
                                   'CLs')) %dopar% {getAccuracy(a)}
  allAccuracy
}

getSingleAccuracyAll <- function(i,Q){
  # Simulate data under DINA model
  gen = DINAsim(Profiles,Q,ss,gs)

  logLikelihood <- t(apply(gen$Y, 1, function(m) colSums(m*log(I1) + (1-m)*log(I2))))
  # Maximum Likelihood Estimation
  MLE <- apply(logLikelihood,1,which.max)
  
  return(sum(MLE==CLs))
}


getSingleAccuracy <- function(i,Q){
  # Simulate data under DINA model
  gen = DINAsim(Profiles,Q,ss,gs)
  
  logLikelihood <- t(apply(gen$Y, 1, function(m) colSums(m*log(I1) + (1-m)*log(I2))))
  # Maximum Likelihood Estimation
  MLE <- apply(logLikelihood,1,which.max)
  result <- sapply(1:8,function(i){
    indexMLE <- MLE==i
    indexCLs <- CLs==i
    temp1 <- CLs[indexMLE] == MLE[indexMLE]
    temp2 <- CLs[indexCLs] != MLE[indexCLs] 
    TP <- sum(temp1)
    FN <- length(temp1)-TP
    FP <- sum(temp2)
    return(c(TP=TP,FN=FN,FP=FP))
  })
  return(as.vector(result))
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
  
  avg <- rowMeans(sapply(1:100, getSingleAccuracy, Q=Q))
  return(c(avg, myIndex))
}  



getCounts <- function(myIndex){
  Q <- Q.complete[myIndex,]
  
  # Simulate data under DINA model
  gen = DINAsim(Profiles,Q,ss,gs)
  
  # Get ideal reponse 
  toCompare <- As %*% t(Q)
  correct <- t(matrix(rep(apply(Q,1,function(x) sum(x^2)), nrow(As)), nrow(Q), nrow(As)))
  idealResponse <- apply(toCompare == correct,2,as.integer)
  
  # Calculate logLikelihood
  I1 <- ((1-ss)^t(idealResponse) * (gs)^t(1-idealResponse))
  I2 <- (ss^t(idealResponse) * (1-gs)^t(1-idealResponse))
  logLikelihood <- t(apply(gen$Y, 1, function(m) colSums(m*log(I1) + (1-m)*log(I2))))
  
  # Maximum Likelihood Estimation
  MLE <- apply(logLikelihood,1,which.max)
  result <- sapply(1:8,function(i){
              indexMLE <- MLE==i
              indexCLs <- CLs==i
              temp1 <- CLs[indexMLE] == MLE[indexMLE]
              temp2 <- CLs[indexCLs] != MLE[indexCLs] 
              TP <- sum(temp1)
              FN <- length(temp1)-TP
              FP <- sum(temp2)
              return(c(TP=TP,FN=FN,FP=FP))
  })
}

library(foreach)
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

set.seed(1)
theta = 0.01
J = 9
ss = gs = rep(theta,J)
myResult <- getAllAccuracy(J,7)
write(myResult, paste0("simPerPatternResult_theta",theta,"_J",J,".txt"))

set.seed(1)
theta = 0.05
J = 9
ss = gs = rep(theta,J)
myResult <- getAllAccuracy(J,7)
write(myResult, paste0("simPerPatternResult_theta",theta,"_J",J,".txt"))

set.seed(1)
theta = 0.1
J = 9
ss = gs = rep(theta,J)
myResult <- getAllAccuracy(J,7)
write(myResult, paste0("simPerPatternResult_theta",theta,"_J",J,".txt"))

set.seed(1)
theta = 0.2
J = 9
ss = gs = rep(theta,J)
myResult <- getAllAccuracy(J,7)
write(myResult, paste0("simPerPatternResult_theta",theta,"_J",J,".txt"))

set.seed(1)
theta = 0.3
J = 9
ss = gs = rep(theta,J)
myResult <- getAllAccuracy(J,7)
write(myResult, paste0("simPerPatternResult_theta",theta,"_J",J,".txt"))

start_time <- Sys.time()
myResult <- getAllAccuracy(J,7)
end_time <- Sys.time()
end_time - start_time
a <- which.max(myResult[1,])
a
myResult[,a]
write(myResult, paste0("simPerPatternResult_theta",theta,"_J",J,".txt"))

# Read simPerPatternResult
require(data.table)
myResult <- scan("simPerPatternResult_theta0.2_J5.txt")
myResult <- matrix(myResult, 29, length(myResult)/29)  # 29 = 8(categories)*3+5(questions)
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

ggplot(myFinal)+geom_point(data = myFinal[,.(configuration, fscore1)])+geom_point(data = myFinal[,.(configuration, fscore2)], colour='red')

require(ggplot2)
library(reshape2)
d <- melt(myFinal, id.vars="configuration")

# Everything on the same plot
ggplot(d, aes(configuration,value, col=variable)) + 
  geom_point() 

ggplot(d[d[,variable %like% "score"]], aes(configuration,value, col = variable)) + 
  geom_point() + 
  geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "score")]],shape = 18) +
  geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "score")]],shape = 23)


g <- ggplot(d[d[,variable %like% "score"]], aes(configuration,value, shape = variable)) + 
      geom_point() + scale_shape_manual(values=1:9) + 
      geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "score")]],col = "red") +
      geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "score")]],col = "blue")
ggsave(g, file = "theta02J5Fscore.svg")

g <- ggplot(d[d[,variable %like% "recision"]], aes(configuration,value, shape = variable)) + 
      geom_point() + scale_shape_manual(values=1:9) + 
      geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "recision")]],col = "red") +
      geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "recision")]],col = "blue")
ggsave(g, file = "theta02J5Precision.svg")

g <- ggplot(d[d[,variable %like% "total"]], aes(configuration,value, shape = variable)) + 
      geom_point() + scale_shape_manual(values=1:9) + 
      geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "total")]],col = "red") +
      geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "total")]],col = "blue")

ggplot(d[d[,variable %like% "total"]], aes(configuration,value, col = variable)) + 
  geom_point() +  
  geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "total")]],shape = 3) +
  geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "total")]],shape = 7)

ggsave(g, file = "theta02J5Total.svg")


library(tikzDevice)
library(ggplot2)
#For some reason, Rstudio needs to know the time zone...
options(tz="CA")

#Create a .tex file that will contain your plot as vectors
#You need to set the size of your plot here, if you do it in LaTeX, font consistency with the rest of the document will be lost
tikz(file = "theta02J5Total.tex", width = 4, height = 3)
#Simple plot of the dummy data using LaTeX elements
g <- ggplot(d[d[,variable %like% "total"]], aes(configuration,value, shape = variable)) + 
  geom_point() + scale_shape_manual(values=1:9) + 
  geom_point(data = d[d[, (configuration %in% myIndex1) & (variable %like% "total")]],col = "red") +
  geom_point(data = d[d[, (configuration %in% myIndex2) & (variable %like% "total")]],col = "blue")
#This line is only necessary if you want to preview the plot right after compiling
print(g)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()



# myResult <- t(read.table("simResult_theta0.2_J5.txt"))
myResult <- scan("simResult_theta0.2_J6.txt")
myResult <- matrix(myResult, 7, length(myResult)/7)
allIndex <- getIndex(5,7)
myIndex1 <- which(apply(allIndex,2,function(m) all(c(1,2,3) %in% m)))
myIndex2 <- which(apply(allIndex,2,function(m) {all(c(1,2,3) %in% m) & any(c(4,5,6,7) %in% m)}))
myResult[,myIndex2]
myData <- data.frame(configuration=c(1:ncol(myResult)), accuracy=myResult[1,]/2)
ggplot(myData, aes(configuration, accuracy))+ geom_point()+ geom_point(data = myData[myIndex1,], colour='red')+ geom_point(data = myData[myIndex2,], colour='blue')

# devtools::install_github('cdeterman/gpuR', ref = 'develop')
# devtools::install_github('gpuRcore/gpuRcuda', ref = 'develop')
# Sys.setenv(OPENCL_INC = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.1\include\CL")
# Sys.setenv(OPENCL_LIB64 = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.1\lib\x64")

setContext(id = 3L) #Using Nvidia

set.seed(123)
ORDER = 1024

A = matrix(rnorm(ORDER^2), nrow=ORDER)
B = matrix(rnorm(ORDER^2), nrow=ORDER)
gpuA = gpuMatrix(A, type="double")
gpuB = gpuMatrix(B, type="double")

C = A %*% B
gpuC = gpuA %*% gpuB

all(C == gpuC[])

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
foreach(i=1:3) %dopar% sqrt(i)


x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
ptime <- system.time({
   r <- foreach(icount(trials), .combine=cbind) %dopar% {
     ind <- sample(100, 100, replace=TRUE)
     result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
     coefficients(result1)
     }
   })

stime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %do% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
  }
})


a <- matrix(1:20, 4,5)
library(iterators)
foreach(b = split(a,row(a)), .combine = c) %do% sum(b)
