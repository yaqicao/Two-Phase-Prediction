require(osDesign)
require(plyr)
require(MESS)
require(mvtnorm)
setwd("~/Desktop/Postdoc in Penn/Oct,2021/10.12-Github")
source("FUN_DATA_sampling.R")
source("FUN_Estimation.R")
source("FUN_Evaluation.R")

#1. Full-Benchmark
Prev <- "Rare" #"Moder" 
Rho <- "0.5" #"Independent"/"0.3"
n <- 3000
if(Prev == "Rare")
{beta <- log(c(0.03,0.6,1.6,0.6,1.5))}
if(Prev == "Moder")
{beta <- log(c(0.13,0.6,1.6,0.6,1.5))}
numbeta <- length(beta)

epsilon_beta <- abs(beta)/5

K=1000
dat_list <- list(NULL)
ResAuc <- matrix(NA, K, 2)
for(k in 1:K)
{
  dat_list[[k]] <- Dat_gen(n, beta, Rho)
  data.full <- dat_list[[k]]
  fullFit_list <- Dat_fullFit(data.full)
  X<-cbind(rep(1,n),data.full$x1,data.full$x2,data.full$x3,data.full$z)
  auc_list_full <- Dat_auc_full(data.full, fullFit_list,X)
  ResAuc[k,] <- c(auc_list_full$auc,auc_list_full$sq_var_auc)
  if(round(k/10)==k/10)
  {print(k)
    print(c(mean(ResAuc[,1],na.rm=T),mean(ResAuc[,2],na.rm=T)))}
}

# 2. Case-control Design
Prev <- "Rare" #"Moder"/"Rare"
Rho <- "Independent" #/"Independent"/"0.3"/"0.5"
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

epsilon_beta <- abs(beta)/5
epsilon_u <- abs(p0)/5
set.seed(123)
K=1000
dat_cc <- list(NULL)
ResAuc <- matrix(NA, K, 2)
for(k in 1:K)
{
  dat <- Dat_gen(n, beta, Rho)
  dat_cc[[k]] <- Dat_gencc(dat, n, beta, p0, ratio)
  
  dat_mlelist <- Dat_format(dat, dat_cc[[k]], "MLE")
  mleFit_list <- Dat_mleFit(dat_mlelist)
  auc_list_mle <- Dat_auc_mle(dat_mlelist, mleFit_list)

  ResAuc[k,] <- c(auc_list_mle$auc,auc_list_mle$var_auc)
  if(round(k/10)==k/10)
  {print(k)
    print(c(mean(ResAuc[,1],na.rm=T),mean(sqrt(ResAuc[,2]),na.rm=T)))}
}

#3. Balanced Design
Prev <- "Rare" #"Moder"/"Rare"
Rho <- "0.3" #/"Independent"/"0.3"/"0.5"
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
epsilon_beta <- abs(beta)/5
epsilon_u <- abs(p0)/5
set.seed(123)
K=1000
dat_Balanced <- list(NULL); dat_BalancedG <- list(NULL)
ResAuc <- matrix(NA, K, 2)
for(k in 1:K)
{
  dat <- Dat_gen(n, beta, Rho)
  G <- sapply(dat$x2, function(x) if (x<=0.25) x<-1 else if (x>0.25 & x<=0.5) x<-2 else if (x>0.5 & x<=0.75) x<-3
              else x<-4)
  dat <- data.frame(cbind(dat,G))
  dat_BalancedG[[k]] <- dat
  dat_Balanced[[k]] <- Dat_genBalanced(dat, n, beta, numG, p0, ratio)
  
  dat_mlelist <- Dat_formatG(dat, dat_Balanced[[k]], "MLE")
  mleFit_list <- Dat_mleFitG(dat_mlelist)
  auc_list_mle <- Dat_auc_mle(dat_mlelist, mleFit_list)

  ResAuc[k,] <- c(auc_list_mle$auc,auc_list_mle$var_auc)
  if(round(k/10)==k/10)
  {print(k)
    print(c(mean(ResAuc[,1],na.rm=T),mean(sqrt(ResAuc[,2]),na.rm=T)))}
}

#4.R-balanced Design
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
ddat <- data.frame(y,x1,x2,x3,z)
glm.fit <- glm(y~1+x1+x2+x3,data=ddat,family = binomial(link = "logit"))
theta <- glm.fit$coefficients

prob <- exp(X[,-5] %*% theta)/(1+exp(X[,-5] %*% theta))
Quan <- unname(quantile(prob[ddat$y==1],c(0.25,0.5,0.75)))
#Quan
G <- sapply(prob, function(x) if (x<=Quan[1]) x<-1 else if (x>Quan[1] & x<=Quan[2]) x<-2 else if (x>Quan[2] & x<=Quan[3]) x<-3
            else x<-4)
dat <- data.frame(cbind(y,x1,x2,x3,z,G))


# Pre-specified number of subjects in phase II by stratum;
x0 <- xtabs(~dat$y+dat$G)[2,]*ratio
p0 <- round(x0/xtabs(~dat$y+dat$G)[1,],3)

epsilon_beta <- abs(beta)/5
epsilon_u <- abs(p0)/5
set.seed(123)
K=1000
dat_Eb <- list(NULL); dat_EbG <- list(NULL)
ResAuc <- matrix(NA, K, 2)
theta.out <- matrix(NA, K, 4); Quan.out <- matrix(NA, K, 3)
for(k in 1:K)
{
  dat <- Dat_gen(n, beta, Rho)
  glm.fit <- glm(y~1+x1+x2+x3,data=dat,family = binomial(link = "logit"))
  theta <- glm.fit$coefficients
  theta.out[k,] <- theta
  X <- cbind(rep(1,n), dat$x1, dat$x2, dat$x3, dat$z)
  prob <- exp(X[,-5] %*% theta)/(1+exp(X[,-5] %*% theta))
  Quan <- unname(quantile(prob[dat$y==1],c(0.25,0.5,0.75)))
  Quan.out[k,] <- Quan
  G <- sapply(prob, function(x) if (x<=Quan[1]) x<-1 else if (x>Quan[1] & x<=Quan[2]) x<-2 else if (x>Quan[2] & x<=Quan[3]) x<-3
              else x<-4)
  
  dat <- data.frame(cbind(dat,G))
  dat_EbG[[k]] <- dat
  dat_Eb[[k]] <- Dat_genEb(dat, n, beta, numG, p0, ratio)
  
  dat_mlelist <- Dat_formatG(dat, dat_Eb[[k]], "MLE")
  mleFit_list <- Dat_mleFitG(dat_mlelist)
  auc_list_mle <- Dat_auc_mle(dat_mlelist, mleFit_list)
  ResAuc[k,] <- c(auc_list_mle$auc,auc_list_mle$var_auc)
  if(round(k/10)==k/10)
  {print(k)
    print(c(mean(ResAuc[,1],na.rm=T),mean(sqrt(ResAuc[,2]),na.rm=T)))}
}
