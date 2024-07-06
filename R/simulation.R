library(matrixcalc)
library(glmnet)
library(pracma)
library(mvtnorm )
library(MCMCpack)
library(BayesLogit)

source("func.R")
load("image.dat")
load("true_coeff.dat")

wv_new <- dim(wv1)[2] #"wv1" comes from image data
nind = 71
ntrain <- 60
ngroup <- 2
wavelet <- wv_new
cov <- 2
para_int <- (wavelet*cov)
para <- wavelet + para_int
###################### Data #######################
x_full <- matrix(0,nrow=nind,ncol=cov)
for(i in 1:nind){
  for(j in 1:cov){
    x_full[i,j] <- rnorm(1,0,0.1)
  }
}
x_train <- x_full[1:ntrain,]

eta_x_true <- rnorm(cov,0,1)
t1 <- rep(0,para_int)
cl1 <- sample(1:para_int,(para_int-4))
length(cl1)
t1[cl1] <- 0
cl2 <- sample(c(1:para_int)[-cl1],2)
t1[cl2] <- -3
cl3 <- c(1:para_int)[-c(cl1,cl2)]
t1[cl3] <- 3
alpha_true <- c(beta_1,t1)#beta1 denotes true regression coefficient


prior_mu_var <- 0.1
prior_beta_var <- 1
prior_sig_alpha <- 1
prior_sig_beta <- 1

C_full <- wv1
D_int_full <- matrix(0,nrow=nind,ncol=para_int)

for(i in 1:nind){
  a <- NULL
  for(j in 1:cov){
    a <- c(a,x_full[i,j]*C_full[i,])
  }
  D_int_full[i,] <- a
}

D_full <- cbind(C_full,D_int_full) 

D_train <- D_full[1:ntrain,]


beta_true <-  rnorm(1,0,0.001)
beta_true
beta_true_vec <- rep(beta_true,nind)
sigma_true <- 0.3

ep <- rnorm(nind,0,(sigma_true)^2)

y_full <- beta_true_vec + (x_full %*% eta_x_true) + (D_full %*%alpha_true) + ep

y <- y_full[1:ntrain,]

mix_prob_prior <- c(0.5,0.5)

error_sig2 <- sigma_true^2

eta_x <- rep(1,cov)

par1 <- matrix(0,nrow=ntrain,ncol=ngroup)
par1_compli <- matrix(0,nrow=ntrain,ncol=ngroup)
for(i in 1:ntrain){
  for(h in 1:ngroup){
    par1[i,h] <- 1
    par1_compli[i,h] <- 1
  }
}
mix_prob <- c(0.5,0.5)
a1 <- rep(0,round(6*para/10))
a2 <- rnorm(round(para/5),-3,0.01)
a3 <- rnorm(para-length(c(a1,a2)),3,0.01)

alpha <- sample(c(a1,a2,a3),para)

alpha_group <- c(0,-3,3)
mu <- 1 #prior mean of alpha 
sigma2 <- 1 #prior variance of alpha

beta <- 1

prior_etax_var <- 0.3*diag(cov) 
#=======================================

iteration <- 5000
burn <- 2500

count <- 0

mix_prob_sim  <- array(0, c(ngroup,1, (iteration-burn)))
alpha_gr_sim <- array(0, c((ngroup+1), 1, (iteration-burn)))
alpha_sim <- array(0, c(para, 1, (iteration-burn)))
beta_sim <- array(0, c(1, 1, (iteration-burn)))
eta_x_sim <- array(0, c(cov, 1, (iteration-burn)))
pred_sim <- array(0, c( (nind-ntrain), 1, (iteration-burn)))
pred_sim_lasso <- array(0, c( (nind-ntrain), 1, B))

start <- proc.time()[3]
max1 <- function(x){
  which(x==max(x))[1]
}

zero_fn <- function(x){
  rbinom(1,1,(1- x))
}

normalize <- function(x){
  if (max(x) < -750 | max(x) > 700) x <- x + (700 - max(x))
  x
}

q <- 0.4
err_alpha <- 1
err_beta <- 1
###################################################

for(ind in 1:iteration){
  
  beta <- rnorm(1,(sum(y- (x_train %*% eta_x) - (D_train %*%alpha) )*prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)),(error_sig2 * prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)) )
  mm1 <- solve((t(x_train)%*% ((1/error_sig2)*diag(ntrain)) %*% x_train ) + solve(prior_etax_var)) %*% 
    t(x_train)%*% ((1/error_sig2)*diag(ntrain))%*%(y - beta - (D_train %*%alpha) )
  mm2 <- solve((t(x_train)%*% ((1/error_sig2)*diag(ntrain)) %*% x_train ) + solve(prior_etax_var))
  eta_x <- mvrnorm(1,mm1,mm2)  
 
  alpha <- alpha_draw(alpha,beta,eta_x,error_sig2,alpha_group,q,mix_prob)
  
  q <- rbeta(1,(para-length(which(alpha==0))) + 1,length(which(alpha==0)) + 1 )
  
  temp1 <- numeric()
  for(h in 1:(ngroup)){
    temp1[h] <- length(which(alpha==alpha_group[h+1]))
  }
  
  mix_prob <- as.vector(rdirichlet(1,(temp1+mix_prob_prior)))
  
  alpha_group <- clus_draw(alpha,alpha_group,beta,eta_x,error_sig2,mu,sigma2)
  
  mu <-  rnorm(1,solve(1+(length(unique(class_alpha))/sigma2))%*%(1+ sum(alpha_group[unique(class_alpha)]/sigma2)),solve(1+(length(unique(class_alpha))/sigma2)))
  sigma2 <- rinvgamma(1,shape=(length(unique(class_alpha))/2)+prior_sig_alpha,scale= (sum((alpha_group[unique(class_alpha)] - mu)^2)/2)+prior_sig_beta  )  
  error_sig2 <- rinvgamma(1,shape=(ntrain/2)+err_alpha,scale= (sum((y-rep(beta,ntrain)- (x_train %*% eta_x) - (D_train%*%alpha))^2)/2 )+ err_beta)
  
  if (ind > burn) 
  {
    count <- count + 1
    y_pred <- numeric()
    for(i in (ntrain+1):nind){
      y_pred[i-ntrain] <- beta + sum(D_full[i,]*alpha)+ sum(x_full[i,]*eta_x)
    }
    pred_sim [ , , count] <- y_pred
    
    mix_prob_sim[,,count]  <- mix_prob
    beta_sim[,,count] <- beta
    alpha_sim[,,count] <- alpha
    alpha_gr_sim[,,count] <- alpha_group
    eta_x_sim[,,count] <- eta_x
  }
  print(ind)
  print(proc.time()[3])
} #iteration

end <- proc.time()[3]
tt <- end - start
tt

zero <- function(x)
{
  z <- length(which(x==0))  
  if (z > (length(x)/2))
  {
    p <- 0
  }else{
    p <- sum(x)/(length(x) - length(which(x==0)))  
  }
  p
}


alpha_est <- apply(alpha_sim[ , , ], 1,  FUN=function(x) zero(x) )
estimation_err <- sqrt(sum((alpha_true - alpha_est)^2)) /sqrt(sum(alpha_true^2))
  
beta_est <- apply(beta_sim, 1, mean)
eta_x_est <- apply(eta_x_sim, 1, mean)
#Prediction
y_pred <- numeric()
for(i in (ntrain+1):nind){
  y_pred[i-ntrain] <- beta_est + sum(D_full[i,]*alpha_est)+ sum(x_full[i,]*eta_x_est)
}
y_test <- y_full[(ntrain+1):nind]

pred_err <- sqrt(sum((y_test - y_pred)^2)) /sqrt(sum((y_test - mean(y_test))^2))

forecast <- forecast_cred(pred_sim,y_test)
cred <- forecast$cred_pred
covergae_prop <- length(which(cred[,3]<cred[,1] & cred[,1]<cred[,7]))/length(cred[,3])
width <- cred[,7]-cred[,3]
