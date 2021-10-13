require(osDesign)
require(plyr)
require(MESS)
require(mvtnorm)
setwd("~/Desktop/Postdoc in Penn/Oct,2021/10.12-Github")
source("FUN_DATA_sampling.R")
source("FUN_Estimation.R")
source("FUN_Evaluation.R")

#CC Design + R-Balanced (Balanced) Post-stratification analysis method

Prev <- "Moder" #"Moder"/"Rare"
Rho <- "0.5" #/"Independent"/"0.3"/"0.5"
ratio <- 2
n <- 3000
if(Prev == "Rare")
  {beta <- log(c(0.03,0.6,1.6,0.6,1.5))}
if(Prev == "Moder")
  {beta <- log(c(0.13,0.6,1.6,0.6,1.5))}
numbeta <- length(beta)
      
x1 <- rnorm(n)
x2 <- runif(n)
x3 <- rbinom(n, size=1, prob=0.2)
z <- rnorm(n)
X <- cbind(rep(1,n), x1, x2, x3, z)
y <- rbinom(n, size=1, prob=exp(X %*% beta)/(1+exp(X %*% beta)))
dat <- data.frame(cbind(y,x1,x2,x3,z))
      
# Pre-specified number of subjects in phase II by stratum;
x0 <- sum(y)*ratio
p0 <- round(x0/sum(1-y), 3)

numG <- 4; numbeta <- 5; n <- 3000
epsilon_beta <- abs(beta)/5
K <- 1000
ResAuc <- matrix(NA, K, 2)
out.beta.mle <- out.sd.tps.covm <- matrix(NA, 1000, numbeta)
for (k in 1:K)
{
  dat<- Dat_gen(n, beta, Rho)
  dat2 <- Dat_gencc(dat, n, beta, p0, ratio)
 #R-balanced Post-stratification analysis method       
  glm.fit <- glm(y~1+x1+x2+x3,data=dat,family = binomial(link = "logit"))
  theta <- glm.fit$coefficients
        
  X2 <- cbind(rep(1,length(dat2$x1)), dat2$x1, dat2$x2, dat2$x3, dat2$z)
  prob2 <- exp(X2[,-5] %*% theta)/(1+exp(X2[,-5] %*% theta))
  Quan <- unname(quantile(prob2,c(0.25,0.5,0.75)))
  G <- sapply(prob2, function(x) if (x<=Quan[1]) x<-1 else if (x>Quan[1] & x<=Quan[2]) x<-2 else if (x>Quan[2] & x<=Quan[3]) x<-3
                    else x<-4)
  dat2 <- cbind(dat2,G)
        
  X <- cbind(rep(1,n), dat$x1, dat$x2, dat$x3, dat$z)
  prob <- exp(X[,-5] %*% theta)/(1+exp(X[,-5] %*% theta))
        
  G <- sapply(prob, function(x) if (x<=Quan[1]) x<-1 else if (x>Quan[1] & x<=Quan[2]) x<-2 else if (x>Quan[2] & x<=Quan[3]) x<-3
                    else x<-4)
        
  dat <- cbind(dat,G)

  ##Balanced Post-stratification analysis method
  #Quan <- unname(quantile(dat2$x2,c(0.25,0.5,0.75)))
  #G <- sapply(dat$x2, function(x) if (x<=Quan[1]) x<-1 else if (x>Quan[1] & x<=Quan[2]) x<-2 else if (x>Quan[2] & x<=Quan[3]) x<-3
  #            else x<-4)
  #dat <- cbind(dat,G)
  
  #G <- sapply(dat2$x2, function(x) if (x<=Quan[1]) x<-1 else if (x>Quan[1] & x<=Quan[2]) x<-2 else if (x>Quan[2] & x<=Quan[3]) x<-3
  #            else x<-4)
  #dat2 <- cbind(dat2,G)
  #xtabs(~dat2$y+dat2$G)
  
  dat_mlelist <- Dat_formatG(dat, dat2, "MLE")
  mleFit_list <- Dat_mleFitG(dat_mlelist)
  auc_list_mle <- Dat_auc_mle(dat_mlelist, mleFit_list)
        
  out.beta.mle[k, ] <- mleFit_list$beta.mle
  out.sd.tps.covm[k, ] <- sqrt(mleFit_list$beta.var.mle)
  ResAuc[k,] <- c(auc_list_mle$auc,sqrt(auc_list_mle$var_auc))
  if(round(k/10)==k/10)
  {
    print(k)
    print(colMeans(ResAuc,na.rm=T))
    print(colMeans(out.beta.mle,na.rm=T))
    print(colMeans(out.sd.tps.covm,na.rm=T))}}
       
  
#Balanced Design + R-Balanced Post-stratification analysis method      
Prev <- "Rare" #"Moder"/"Rare"
Rho <- "Independent" #/"Independent"/"0.3"/"0.5"
ratio <- 2
set.seed(719)
n <- 3000
if(Prev == "Rare")
{beta <- log(c(0.03,0.6,1.6,0.6,1.5))}
if(Prev == "Moder")
{beta <- log(c(0.13,0.6,1.6,0.6,1.5))}
numbeta <- length(beta)
numG <- 4
x1 <- rnorm(n)
x2 <- runif(n)
x3 <- rbinom(n, size=1, prob=0.2)
z <- rnorm(n)
X <- cbind(rep(1,n), x1, x2, x3, z)
y <- rbinom(n, size=1, prob=exp(X %*% beta)/(1+exp(X %*% beta)))
G <- sapply(x2, function(x) if (x<=0.25) x<-1 else if (x>0.25 & x<=0.5) x<-2 else if (x>0.5 & x<=0.75) x<-3
            else x<-4)
dat <- data.frame(cbind(y,x1,x2,x3,z,G))


# Pre-specified number of subjects in phase II by stratum;
x0 <- xtabs(~dat$y+dat$G)[2,]*ratio
p0 <- round(x0/xtabs(~dat$y+dat$G)[1,],3)       

numP <- 2 
epsilon_beta <- abs(beta)/5
K <- 1000
ResAuc <- matrix(NA, K, 2)
out.beta.mle <- out.sd.tps.covm <- out.sd.tps.cove <- matrix(NA, 1000, numbeta)
for (k in 1:K)
  {numG <- 4
  dat<- Dat_gen(n, beta, Rho)
  G <- sapply(dat$x2, function(x) if (x<=0.25) x<-1 else if (x>0.25 & x<=0.5) x<-2 else if (x>0.5 & x<=0.75) x<-3
              else x<-4)
  dat <- data.frame(cbind(dat,G))
  dat2 <- Dat_genBalanced(dat, n, beta, numG, p0, ratio)
            
  glm.fit <- glm(y~1+x1+x2+x3,data=dat,family = binomial(link = "logit"))
  theta <- glm.fit$coefficients
            
  X2 <- cbind(rep(1,length(dat2$x1)), dat2$x1, dat2$x2, dat2$x3, dat2$z)
  dat2Star <- NULL;Quan <- matrix(NA, numG, numP-1)
  for (ii in 1:numG)
     {XX2 <- X2[which(dat2$G==ii),]
     prob2 <- exp(XX2[,-5] %*% theta)/(1+exp(XX2[,-5] %*% theta))
     Quan[ii,] <- unname(quantile(prob2,c(1:(numP-1))/numP))
     if (numP==2){
         Gstar <- (ii-1)*numP+sapply(prob2, function(x) if (x<=Quan[ii,1]) x<-1 else x<-2)}
         dat2Star <- rbind(dat2Star,cbind(dat2[which(dat2$G==ii),-6],Gstar))}
         #xtabs(~dat2Star$y+dat2Star$Gstar)
            
        X <- cbind(rep(1,n), dat$x1, dat$x2, dat$x3, dat$z)
        datStar <- NULL
        for (ii in 1:numG)
        {XX <- X[which(dat$G==ii),]
        prob <- exp(XX[,-5] %*% theta)/(1+exp(XX[,-5] %*% theta))
        
        if(numP==2){
        Gstar <- (ii-1)*numP+sapply(prob, function(x) if (x<=Quan[ii,1]) x<-1 else x<-2)}
        datStar <- rbind(datStar,cbind(dat[which(dat$G==ii),-6],Gstar))}
        #xtabs(~dat2Star$y+dat2Star$G) 
        #xtabs(~datStar$y+datStar$G)
        #merge(datStar,dat2Star)
        colnames(datStar) <-  c("y","x1","x2","x3","z","G")
        colnames(dat2Star) <-  c("y","x1","x2","x3","z","G")
        numG <- numG*numP
        dat_mlelist <- Dat_formatG(datStar, dat2Star, "MLE")
        data_list <- dat_mlelist
        mleFit_list <- Dat_mleFitG(dat_mlelist)
        auc_list_mle <- Dat_auc_mle(dat_mlelist, mleFit_list)
            
        out.beta.mle[k, ] <- mleFit_list$beta.mle
        out.sd.tps.covm[k, ] <- sqrt(mleFit_list$beta.var.mle)
           
        ResAuc[k,] <- c(auc_list_mle$auc,sqrt(auc_list_mle$var_auc))
        if(round(k/10)==k/10)
         {print(k)
          print(colMeans(ResAuc,na.rm=T))
          print(colMeans(out.beta.mle,na.rm=T))
          print(colMeans(out.sd.tps.covm,na.rm=T))}}
