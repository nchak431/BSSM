#Evaluates the prediction performance of BSSM model
source("func.R")
library(wavethresh)
library(matrixcalc)
library(glmnet)
library(pracma)
library(mvtnorm )
library(MCMCpack)
library(BayesLogit)
load("image.dat") #image data after wavelet decomposition
load("res_cov.dat") # response and covariate data
load("index.dat") #marks the indices corresponding to non-zero wavelet coefficients of the image

dim <- 64
dims=c(64,64)
family = "DaubLeAsymm"
fil.num = 8 
level = 3
wv_new <- dim(wv1)[2] #"wv1" comes from image data
nind = 71
ntrain <- 60
ngroup <- 2
wavelet <- wv_new
cov <- 4
para_int <- (wavelet*cov)
para <- wavelet + para_int
###################### Data #######################
x_full <- as.matrix(res_cov_bl[,3:6]) # comes from covariate data
x_train <- as.matrix(x_full[1:ntrain,])

eta_x_true <- rnorm(cov,0,1)

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

y_full <- as.matrix(res_cov_bl[,8]) # comes from response data

y <- y_full[1:ntrain,]

mix_prob_prior <- c(0.5,0.5)
error_sig2 <- 0.3^2
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
a1 <- rep(0,round(9*para/10))
a2 <- rnorm(round(para/20),-0.5,0.01)
a3 <- rnorm(para-length(c(a1,a2)),0.5,0.01)
alpha <- sample(c(a1,a2,a3),para)

alpha_group <- c(0,-0.5,0.5)
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
inv_img_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int1_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int2_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int3_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int4_sim <- array(0, c(dim, dim, (iteration-burn)))
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

q <- 0.5
err_alpha <- 1
err_beta <- 1
###################################################

for(ind in 1:iteration){
  alpha_img <- rep(0,dim^2)
  alpha_int1 <- rep(0,dim^2)
  alpha_int2 <- rep(0,dim^2)
  alpha_int3 <- rep(0,dim^2)
  alpha_int4 <- rep(0,dim^2)
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
  alpha_img[ind_c1] <- alpha[1:wavelet] #ind_c1 comes from the index data
  alpha_int1[ind_c1] <- alpha[(wavelet+1):(2*wavelet)]
  alpha_int2[ind_c1] <- alpha[((2*wavelet)+1):(3*wavelet)]
  alpha_int3[ind_c1] <- alpha[((3*wavelet)+1):(4*wavelet)]
  alpha_int4[ind_c1] <- alpha[((4*wavelet)+1):para]
  
  inv_img <- reconstr.2D(alpha_img,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int1 <- reconstr.2D(alpha_int1,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int2 <- reconstr.2D(alpha_int2,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int3 <- reconstr.2D(alpha_int3,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int4 <- reconstr.2D(alpha_int4,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)#matrix.heatmap(inv_img)
  
  if (ind > burn)
  {
    count <- count + 1 
    mix_prob_sim[,,count]  <- mix_prob
    beta_sim[,,count] <- beta
    alpha_sim[,,count] <- alpha
    inv_img_sim [,,count] <- inv_img
    inv_int1_sim [,,count] <- inv_int1
    inv_int2_sim [,,count] <- inv_int2
    inv_int3_sim [,,count] <- inv_int3
    inv_int4_sim [,,count] <- inv_int4
    
    alpha_gr_sim[,,count] <- alpha_group
    eta_x_sim[,,count] <- eta_x
  }
  print(ind)
  print(table(alpha))
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

### Sig voxels

cred_alpha <- cred_voxel(inv_img_sim)$cred
cred_int1 <- cred_voxel(inv_int1_sim)$cred
cred_int2 <- cred_voxel(inv_int2_sim)$cred
cred_int3 <- cred_voxel(inv_int3_sim)$cred
cred_int4 <- cred_voxel(inv_int4_sim)$cred

################### BCM
library(horseshoe)

design_mat <- cbind(rep(1,ntrain),x_train,D_train)
res <- horseshoe(y, design_mat, nmc=4500,method.tau = "truncatedCauchy", method.sigma = "Jeffreys")
#res$TauHat #posterior mean of tau
coef <- res$BetaHat
samples <- res$BetaSamples
dim(samples)


inv_img_bcm_sim <- array(0, c(dim, dim, 4500))

for(i in 1:4500){
  alpha_img <- rep(0,dim^2)
  alpha_img[ind_c1] <- samples[,i][(cov+2):(wavelet+1+cov)]
  inv_img <- reconstr.2D(alpha_img,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_img_bcm_sim [,,i] <- inv_img
}
cred_alpha_bcm <- cred_voxel(inv_img_bcm_sim)$cred

##################

save(cred_alpha, file=paste("/home/nchakraborty/blue/clus/model2_int",
                            "/real_data/batch2/res_ad1/124_cred_alpha_test-", index1, ".dat", sep=''))

save(cred_int1, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/124_cred_int1_test-", index1, ".dat", sep=''))
save(cred_int2, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/124_cred_int2_test-", index1, ".dat", sep=''))
save(cred_int3, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/124_cred_int3_test-", index1, ".dat", sep=''))
save(cred_int4, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/124_cred_int4_test-", index1, ".dat", sep=''))

save(cred_alpha_bcm, file=paste("/home/nchakraborty/blue/clus/model2_int",
                                "/real_data/batch2/res_ad1/124_cred_alpha_bcm_test-", index1, ".dat", sep=''))
                                
###################################

################################################
load("wv_sl125_part1_ad.dat")
load("res_cov_ad_bl.dat")
load("ind_c1_125_part1_ad.dat")

dim <- 64
dims=c(64,64)
family = "DaubLeAsymm"
fil.num = 8 #Can use 4/6/8/10
level = 3
wv_new <- dim(wv1)[2]
nind = 71
ntrain <- 71
ngroup <- 2
wavelet <- wv_new
cov <- 4
para_int <- (wavelet*cov)
para <- wavelet + para_int
###################### Data #######################
x_full <- as.matrix(res_cov_bl[,3:6])
x_train <- as.matrix(x_full[1:ntrain,])

eta_x_true <- rnorm(cov,0,1)

prior_mu_var <- 0.1
prior_beta_var <- 1
prior_sig_alpha <- 1
prior_sig_beta <- 1

# ================ Data =============
# C_full <- matrix(0,nrow=nind,ncol=wavelet)
#
# for(i in 1:nind){
#   for(k in 1:wavelet){
#     C_full[i,k] <- rnorm(1,0.001,0.1)
#   }
# }
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



y_full <- as.matrix(res_cov_bl[,8])

y <- y_full[1:ntrain,]

mix_prob_prior <- c(0.5,0.5)
############## Initialize #####################
#spike_ind: binary indicator variable
error_sig2 <- 0.3^2

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

# x <- cbind(D_train,x_train)
# cv_model <- cv.glmnet(x, y, alpha = 1)
# best_lambda <- cv_model$lambda.min
# best_lambda
# best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
# coeff_lasso <- as.vector(coef(best_model))
# alpha <- coeff_lasso[2:(para+1)]
# df <- kmeans(alpha,3)
# alpha_group <- as.vector(df$centers)

#a1 <- rnorm(round(9*para/10),0,0.001)
a1 <- rep(0,round(9*para/10))
a2 <- rnorm(round(para/20),-0.5,0.01)
a3 <- rnorm(para-length(c(a1,a2)),0.5,0.01)

alpha <- sample(c(a1,a2,a3),para)

alpha_group <- c(0,-0.5,0.5)
mu <- 1 #prior mean of alpha
sigma2 <- 1 #prior variance of alpha

#
# reg <- lm(y~D_train)
# out <- reg$coefficients
# df <- kmeans(out,2)
#
# alpha <- out[-1]
# alpha_group <- as.vector(df$centers)
# mu <- alpha_group #prior mean of alpha
# sigma2 <- rep(0.1,ngroup) #prior variance of alpha

beta <- 1

prior_etax_var <- 0.3*diag(cov)
B=500
#=======================================

iteration <- 1500
burn <- 700

count <- 0

mix_prob_sim  <- array(0, c(ngroup,1, (iteration-burn)))
alpha_gr_sim <- array(0, c((ngroup+1), 1, (iteration-burn)))
alpha_sim <- array(0, c(para, 1, (iteration-burn)))
inv_img_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int1_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int2_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int3_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int4_sim <- array(0, c(dim, dim, (iteration-burn)))
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

q <- 0.5# This is important
err_alpha <- 1
err_beta <- 1
###################################################

for(ind in 1:iteration){
  alpha_img <- rep(0,dim^2)
  alpha_int1 <- rep(0,dim^2)
  alpha_int2 <- rep(0,dim^2)
  alpha_int3 <- rep(0,dim^2)
  alpha_int4 <- rep(0,dim^2)
  beta <- rnorm(1,(sum(y- (x_train %*% eta_x) - (D_train %*%alpha) )*prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)),(error_sig2 * prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)) )
  # beta <- rnorm(1,(sum(y- (x_train * eta_x) - (D_train %*%alpha) )*prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)),(error_sig2 * prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)) )
  mm1 <- solve((t(x_train)%*% ((1/error_sig2)*diag(ntrain)) %*% x_train ) + solve(prior_etax_var)) %*%
    t(x_train)%*% ((1/error_sig2)*diag(ntrain))%*%(y - beta - (D_train %*%alpha) )
  mm2 <- solve((t(x_train)%*% ((1/error_sig2)*diag(ntrain)) %*% x_train ) + solve(prior_etax_var))
  eta_x <- mvrnorm(1,mm1,mm2)
  # for(h in 1:ngroup){
  #   sum1 <- NULL
  #
  #   for(k in 1:para){
  #
  #     # tmp1 <- numeric()
  #     # for(i in 1:ntrain){
  #     #   tmp1[i] <- sum(alpha*D_train[i,]) - (alpha[k]*D_train[i,k])
  #     # }
  #     tmp1 <- D_train%*%alpha - alpha[k]*D_train[,k]
  #
  #     sum1[k] <- mix_prob[h]*exp( (-1/(2*error_sig2)) *
  #                                   (((alpha_group[h]^2)*sum(D_train[,k]^2))- (2*alpha_group[h]*sum(D_train[,k]*(y- beta-(x_train %*% eta_x) -tmp1))) ) )
  #   }
  #   gr[,h] <- sum1
  # }
  
  for(k in 1:para){
    gr <- numeric()
    # c1 <- matrix(rep(D_train%*%alpha,para),nrow=ntrain)-hadamard.prod(matrix(rep(alpha,ntrain),nrow=ntrain,byrow = T),D_train)
    c2 <- sum(D_train[,k]^2)
    c3 <- sum(D_train[,k]*(y- beta-(x_train %*% eta_x)-((D_train%*%alpha)-(D_train[,k]*alpha[k]))))
    #c3 <- alpha_group[h+1]*apply(hadamard.prod(D_train,matrix(rep(y- beta-(x_train * eta_x),para),nrow=ntrain) -c1),2,sum  )
    gr <- (-1/(2*error_sig2)) *(((alpha_group[2:(ngroup+1)]^2)*c2)-(2*alpha_group[2:(ngroup+1)]*c3))
    gr <- c(0,gr)
    if (max(gr) < -750 | max(gr) > 700) gr <- gr + (700 - max(gr))  
    #gr1 <- t(as.matrix(apply(gr, 1, FUN=function(x) normalize(x) )))
    
    mix1 <- c((1-q),  q*mix_prob)
    
    prob <- hadamard.prod(exp(gr),mix1)
    
    class_alpha <- which(prob==max(prob))[1]
    
    alpha[k] <- alpha_group[class_alpha]
    #print(k)
  }
  
  
  q <- rbeta(1,(para-length(which(alpha==0))) + 1,length(which(alpha==0)) + 1 )
  
  # for(k in 1:para){
  #    class <- which(gr[k,]==max(gr[k,]))
  #    alpha[k] <- alpha_group[class]
  #  }
  
  temp1 <- numeric()
  for(h in 1:(ngroup)){
    temp1[h] <- length(which(alpha==alpha_group[h+1]))
  }
  
  mix_prob <- as.vector(rdirichlet(1,(temp1+mix_prob_prior)))
  
  for(h in 1:(ngroup)){
    if(temp1[h]>0){
      for(i in 1:ntrain){
        index <-  which(alpha==alpha_group[h+1])
        index_compli <-  which(alpha!=alpha_group[h+1])
        par1[i,h] <- sum(D_train[i,index])
        par1_compli[i,h] <- sum(alpha[index_compli]*D_train[i,index_compli])
      }
      
      a1 <- sum((par1[,h])^2)
      a2 <- sum(par1[,h]*(y-rep(beta,ntrain)- (x_train %*% eta_x)-par1_compli[,h] ))
      #a2 <- sum(par1[,h]*(y-rep(beta,ntrain)- (x_train * eta_x)-par1_compli[,h] ))
      mm <- ((a2/error_sig2) + (mu/sigma2))/((a1/error_sig2) + (1/sigma2))
      ss <- 1/((a1/error_sig2) + (1/sigma2))
      alpha_group[h+1] <- rnorm(1,mm,ss)
    }
    else{
      alpha_group[h+1] <- rnorm(1,mu, sigma2)
      
    }
  }
  mu <-  rnorm(1,solve(1+(length(unique(class_alpha))/sigma2))%*%(1+ sum(alpha_group[unique(class_alpha)]/sigma2)),solve(1+(length(unique(class_alpha))/sigma2)))
  sigma2 <- rinvgamma(1,shape=(length(unique(class_alpha))/2)+prior_sig_alpha,scale= (sum((alpha_group[unique(class_alpha)] - mu)^2)/2)+prior_sig_beta  )
  
  error_sig2 <- rinvgamma(1,shape=(ntrain/2)+err_alpha,scale= (sum((y-rep(beta,ntrain)- (x_train %*% eta_x) - (D_train%*%alpha))^2)/2 )+ err_beta)
  #error_sig2 <- rinvgamma(1,shape=(ntrain/2)+err_alpha,scale= (sum((y-rep(beta,ntrain)- (x_train * eta_x) - (D_train%*%alpha))^2)/2 )+ err_beta)
  alpha_img[ind_c1] <- alpha[1:wavelet]
  alpha_int1[ind_c1] <- alpha[(wavelet+1):(2*wavelet)]
  alpha_int2[ind_c1] <- alpha[((2*wavelet)+1):(3*wavelet)]
  alpha_int3[ind_c1] <- alpha[((3*wavelet)+1):(4*wavelet)]
  alpha_int4[ind_c1] <- alpha[((4*wavelet)+1):para]
  
  inv_img <- reconstr.2D(alpha_img,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int1 <- reconstr.2D(alpha_int1,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int2 <- reconstr.2D(alpha_int2,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int3 <- reconstr.2D(alpha_int3,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int4 <- reconstr.2D(alpha_int4,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)#matrix.heatmap(inv_img)
  #matrix.heatmap(inv_img)
  
  if (ind > burn)
  {
    count <- count + 1
  
    mix_prob_sim[,,count]  <- mix_prob
    beta_sim[,,count] <- beta
    alpha_sim[,,count] <- alpha
    inv_img_sim [,,count] <- inv_img
    inv_int1_sim [,,count] <- inv_int1
    inv_int2_sim [,,count] <- inv_int2
    inv_int3_sim [,,count] <- inv_int3
    inv_int4_sim [,,count] <- inv_int4
    
    alpha_gr_sim[,,count] <- alpha_group
    eta_x_sim[,,count] <- eta_x
  }
  print(ind)
  print(table(alpha))
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

### Sig voxels

cred_alpha <- cred_voxel(inv_img_sim)$cred
cred_int1 <- cred_voxel(inv_int1_sim)$cred
cred_int2 <- cred_voxel(inv_int2_sim)$cred
cred_int3 <- cred_voxel(inv_int3_sim)$cred
cred_int4 <- cred_voxel(inv_int4_sim)$cred

################### BCM
library(horseshoe)

design_mat <- cbind(rep(1,ntrain),x_train,D_train)
res <- horseshoe(y, design_mat, nmc=4500,method.tau = "truncatedCauchy", method.sigma = "Jeffreys")
#res$TauHat #posterior mean of tau
coef <- res$BetaHat
samples <- res$BetaSamples
dim(samples)


inv_img_bcm_sim <- array(0, c(dim, dim, 4500))

for(i in 1:4500){
  alpha_img <- rep(0,dim^2)
  alpha_img[ind_c1] <- samples[,i][(cov+2):(wavelet+1+cov)]
  inv_img <- reconstr.2D(alpha_img,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_img_bcm_sim [,,i] <- inv_img
}
cred_alpha_bcm <- cred_voxel(inv_img_bcm_sim)$cred
##################

save(cred_alpha, file=paste("/home/nchakraborty/blue/clus/model2_int",
                            "/real_data/batch2/res_ad1/125_cred_alpha_test-", index1, ".dat", sep=''))

save(cred_int1, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/125_cred_int1_test-", index1, ".dat", sep=''))
save(cred_int2, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/125_cred_int2_test-", index1, ".dat", sep=''))
save(cred_int3, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/125_cred_int3_test-", index1, ".dat", sep=''))
save(cred_int4, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/125_cred_int4_test-", index1, ".dat", sep=''))

save(cred_alpha_bcm, file=paste("/home/nchakraborty/blue/clus/model2_int",
                                "/real_data/batch2/res_ad1/125_cred_alpha_bcm_test-", index1, ".dat", sep=''))      

###############################

load("wv_sl126_part1_ad.dat")
load("res_cov_ad_bl.dat")
load("ind_c1_126_part1_ad.dat")

dim <- 64
dims=c(64,64)
family = "DaubLeAsymm"
fil.num = 8 #Can use 4/6/8/10
level = 3
wv_new <- dim(wv1)[2]
nind = 71
ntrain <- 71
ngroup <- 2
wavelet <- wv_new
cov <- 4
para_int <- (wavelet*cov)
para <- wavelet + para_int
###################### Data #######################
x_full <- as.matrix(res_cov_bl[,3:6])
x_train <- as.matrix(x_full[1:ntrain,])

eta_x_true <- rnorm(cov,0,1)

prior_mu_var <- 0.1
prior_beta_var <- 1
prior_sig_alpha <- 1
prior_sig_beta <- 1

# ================ Data =============
# C_full <- matrix(0,nrow=nind,ncol=wavelet)
#
# for(i in 1:nind){
#   for(k in 1:wavelet){
#     C_full[i,k] <- rnorm(1,0.001,0.1)
#   }
# }
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



y_full <- as.matrix(res_cov_bl[,8])

y <- y_full[1:ntrain,]

mix_prob_prior <- c(0.5,0.5)
############## Initialize #####################
#spike_ind: binary indicator variable
error_sig2 <- 0.3^2

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

# x <- cbind(D_train,x_train)
# cv_model <- cv.glmnet(x, y, alpha = 1)
# best_lambda <- cv_model$lambda.min
# best_lambda
# best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
# coeff_lasso <- as.vector(coef(best_model))
# alpha <- coeff_lasso[2:(para+1)]
# df <- kmeans(alpha,3)
# alpha_group <- as.vector(df$centers)

#a1 <- rnorm(round(9*para/10),0,0.001)
a1 <- rep(0,round(9*para/10))
a2 <- rnorm(round(para/20),-0.5,0.01)
a3 <- rnorm(para-length(c(a1,a2)),0.5,0.01)

alpha <- sample(c(a1,a2,a3),para)

alpha_group <- c(0,-0.5,0.5)
mu <- 1 #prior mean of alpha
sigma2 <- 1 #prior variance of alpha

#
# reg <- lm(y~D_train)
# out <- reg$coefficients
# df <- kmeans(out,2)
#
# alpha <- out[-1]
# alpha_group <- as.vector(df$centers)
# mu <- alpha_group #prior mean of alpha
# sigma2 <- rep(0.1,ngroup) #prior variance of alpha

beta <- 1

prior_etax_var <- 0.3*diag(cov)
B=500
#=======================================

iteration <- 1500
burn <- 700

count <- 0

mix_prob_sim  <- array(0, c(ngroup,1, (iteration-burn)))
alpha_gr_sim <- array(0, c((ngroup+1), 1, (iteration-burn)))
alpha_sim <- array(0, c(para, 1, (iteration-burn)))
inv_img_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int1_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int2_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int3_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int4_sim <- array(0, c(dim, dim, (iteration-burn)))
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

q <- 0.5# This is important
err_alpha <- 1
err_beta <- 1
###################################################

for(ind in 1:iteration){
  alpha_img <- rep(0,dim^2)
  alpha_int1 <- rep(0,dim^2)
  alpha_int2 <- rep(0,dim^2)
  alpha_int3 <- rep(0,dim^2)
  alpha_int4 <- rep(0,dim^2)
  beta <- rnorm(1,(sum(y- (x_train %*% eta_x) - (D_train %*%alpha) )*prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)),(error_sig2 * prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)) )
  # beta <- rnorm(1,(sum(y- (x_train * eta_x) - (D_train %*%alpha) )*prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)),(error_sig2 * prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)) )
  mm1 <- solve((t(x_train)%*% ((1/error_sig2)*diag(ntrain)) %*% x_train ) + solve(prior_etax_var)) %*%
    t(x_train)%*% ((1/error_sig2)*diag(ntrain))%*%(y - beta - (D_train %*%alpha) )
  mm2 <- solve((t(x_train)%*% ((1/error_sig2)*diag(ntrain)) %*% x_train ) + solve(prior_etax_var))
  eta_x <- mvrnorm(1,mm1,mm2)
  # for(h in 1:ngroup){
  #   sum1 <- NULL
  #
  #   for(k in 1:para){
  #
  #     # tmp1 <- numeric()
  #     # for(i in 1:ntrain){
  #     #   tmp1[i] <- sum(alpha*D_train[i,]) - (alpha[k]*D_train[i,k])
  #     # }
  #     tmp1 <- D_train%*%alpha - alpha[k]*D_train[,k]
  #
  #     sum1[k] <- mix_prob[h]*exp( (-1/(2*error_sig2)) *
  #                                   (((alpha_group[h]^2)*sum(D_train[,k]^2))- (2*alpha_group[h]*sum(D_train[,k]*(y- beta-(x_train %*% eta_x) -tmp1))) ) )
  #   }
  #   gr[,h] <- sum1
  # }
  
  for(k in 1:para){
    gr <- numeric()
    # c1 <- matrix(rep(D_train%*%alpha,para),nrow=ntrain)-hadamard.prod(matrix(rep(alpha,ntrain),nrow=ntrain,byrow = T),D_train)
    c2 <- sum(D_train[,k]^2)
    c3 <- sum(D_train[,k]*(y- beta-(x_train %*% eta_x)-((D_train%*%alpha)-(D_train[,k]*alpha[k]))))
    #c3 <- alpha_group[h+1]*apply(hadamard.prod(D_train,matrix(rep(y- beta-(x_train * eta_x),para),nrow=ntrain) -c1),2,sum  )
    gr <- (-1/(2*error_sig2)) *(((alpha_group[2:(ngroup+1)]^2)*c2)-(2*alpha_group[2:(ngroup+1)]*c3))
    gr <- c(0,gr)
    if (max(gr) < -750 | max(gr) > 700) gr <- gr + (700 - max(gr))  
    #gr1 <- t(as.matrix(apply(gr, 1, FUN=function(x) normalize(x) )))
    
    mix1 <- c((1-q),  q*mix_prob)
    
    prob <- hadamard.prod(exp(gr),mix1)
    
    class_alpha <- which(prob==max(prob))[1]
    
    alpha[k] <- alpha_group[class_alpha]
    #print(k)
  }
  
  
  q <- rbeta(1,(para-length(which(alpha==0))) + 1,length(which(alpha==0)) + 1 )
  
  # for(k in 1:para){
  #    class <- which(gr[k,]==max(gr[k,]))
  #    alpha[k] <- alpha_group[class]
  #  }
  
  temp1 <- numeric()
  for(h in 1:(ngroup)){
    temp1[h] <- length(which(alpha==alpha_group[h+1]))
  }
  
  mix_prob <- as.vector(rdirichlet(1,(temp1+mix_prob_prior)))
  
  for(h in 1:(ngroup)){
    if(temp1[h]>0){
      for(i in 1:ntrain){
        index <-  which(alpha==alpha_group[h+1])
        index_compli <-  which(alpha!=alpha_group[h+1])
        par1[i,h] <- sum(D_train[i,index])
        par1_compli[i,h] <- sum(alpha[index_compli]*D_train[i,index_compli])
      }
      
      a1 <- sum((par1[,h])^2)
      a2 <- sum(par1[,h]*(y-rep(beta,ntrain)- (x_train %*% eta_x)-par1_compli[,h] ))
      #a2 <- sum(par1[,h]*(y-rep(beta,ntrain)- (x_train * eta_x)-par1_compli[,h] ))
      mm <- ((a2/error_sig2) + (mu/sigma2))/((a1/error_sig2) + (1/sigma2))
      ss <- 1/((a1/error_sig2) + (1/sigma2))
      alpha_group[h+1] <- rnorm(1,mm,ss)
    }
    else{
      alpha_group[h+1] <- rnorm(1,mu, sigma2)
      
    }
  }
  mu <-  rnorm(1,solve(1+(length(unique(class_alpha))/sigma2))%*%(1+ sum(alpha_group[unique(class_alpha)]/sigma2)),solve(1+(length(unique(class_alpha))/sigma2)))
  sigma2 <- rinvgamma(1,shape=(length(unique(class_alpha))/2)+prior_sig_alpha,scale= (sum((alpha_group[unique(class_alpha)] - mu)^2)/2)+prior_sig_beta  )
  
  error_sig2 <- rinvgamma(1,shape=(ntrain/2)+err_alpha,scale= (sum((y-rep(beta,ntrain)- (x_train %*% eta_x) - (D_train%*%alpha))^2)/2 )+ err_beta)
  #error_sig2 <- rinvgamma(1,shape=(ntrain/2)+err_alpha,scale= (sum((y-rep(beta,ntrain)- (x_train * eta_x) - (D_train%*%alpha))^2)/2 )+ err_beta)
  alpha_img[ind_c1] <- alpha[1:wavelet]
  alpha_int1[ind_c1] <- alpha[(wavelet+1):(2*wavelet)]
  alpha_int2[ind_c1] <- alpha[((2*wavelet)+1):(3*wavelet)]
  alpha_int3[ind_c1] <- alpha[((3*wavelet)+1):(4*wavelet)]
  alpha_int4[ind_c1] <- alpha[((4*wavelet)+1):para]
  
  inv_img <- reconstr.2D(alpha_img,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int1 <- reconstr.2D(alpha_int1,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int2 <- reconstr.2D(alpha_int2,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int3 <- reconstr.2D(alpha_int3,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int4 <- reconstr.2D(alpha_int4,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)#matrix.heatmap(inv_img)
  #matrix.heatmap(inv_img)
  
  if (ind > burn)
  {
    count <- count + 1
    
    mix_prob_sim[,,count]  <- mix_prob
    beta_sim[,,count] <- beta
    alpha_sim[,,count] <- alpha
    inv_img_sim [,,count] <- inv_img
    inv_int1_sim [,,count] <- inv_int1
    inv_int2_sim [,,count] <- inv_int2
    inv_int3_sim [,,count] <- inv_int3
    inv_int4_sim [,,count] <- inv_int4
    
    alpha_gr_sim[,,count] <- alpha_group
    eta_x_sim[,,count] <- eta_x
  }
  print(ind)
  print(table(alpha))
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

### Sig voxels

cred_alpha <- cred_voxel(inv_img_sim)$cred
cred_int1 <- cred_voxel(inv_int1_sim)$cred
cred_int2 <- cred_voxel(inv_int2_sim)$cred
cred_int3 <- cred_voxel(inv_int3_sim)$cred
cred_int4 <- cred_voxel(inv_int4_sim)$cred

################### BCM
library(horseshoe)

design_mat <- cbind(rep(1,ntrain),x_train,D_train)
res <- horseshoe(y, design_mat, nmc=4500,method.tau = "truncatedCauchy", method.sigma = "Jeffreys")
#res$TauHat #posterior mean of tau
coef <- res$BetaHat
samples <- res$BetaSamples
dim(samples)


inv_img_bcm_sim <- array(0, c(dim, dim, 4500))

for(i in 1:4500){
  alpha_img <- rep(0,dim^2)
  alpha_img[ind_c1] <- samples[,i][(cov+2):(wavelet+1+cov)]
  inv_img <- reconstr.2D(alpha_img,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_img_bcm_sim [,,i] <- inv_img
}
cred_alpha_bcm <- cred_voxel(inv_img_bcm_sim)$cred
##################

save(cred_alpha, file=paste("/home/nchakraborty/blue/clus/model2_int",
                            "/real_data/batch2/res_ad1/126_cred_alpha_test-", index1, ".dat", sep=''))

save(cred_int1, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/126_cred_int1_test-", index1, ".dat", sep=''))
save(cred_int2, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/126_cred_int2_test-", index1, ".dat", sep=''))
save(cred_int3, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/126_cred_int3_test-", index1, ".dat", sep=''))
save(cred_int4, file=paste("/home/nchakraborty/blue/clus/model2_int",
                           "/real_data/batch2/res_ad1/126_cred_int4_test-", index1, ".dat", sep=''))

save(cred_alpha_bcm, file=paste("/home/nchakraborty/blue/clus/model2_int",
                                "/real_data/batch2/res_ad1/126_cred_alpha_bcm_test-", index1, ".dat", sep=''))                                

#############################################

load("wv_sl127_part1_ad.dat")
load("res_cov_ad_bl.dat")
load("ind_c1_127_part1_ad.dat")

dim <- 64
dims=c(64,64)
family = "DaubLeAsymm"
fil.num = 8 #Can use 4/6/8/10
level = 3
wv_new <- dim(wv1)[2]
nind = 71
ntrain <- 71
ngroup <- 2
wavelet <- wv_new
cov <- 4
para_int <- (wavelet*cov)
para <- wavelet + para_int
###################### Data #######################
x_full <- as.matrix(res_cov_bl[,3:6])
x_train <- as.matrix(x_full[1:ntrain,])

eta_x_true <- rnorm(cov,0,1)

prior_mu_var <- 0.1
prior_beta_var <- 1
prior_sig_alpha <- 1
prior_sig_beta <- 1

# ================ Data =============
# C_full <- matrix(0,nrow=nind,ncol=wavelet)
#
# for(i in 1:nind){
#   for(k in 1:wavelet){
#     C_full[i,k] <- rnorm(1,0.001,0.1)
#   }
# }
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



y_full <- as.matrix(res_cov_bl[,8])

y <- y_full[1:ntrain,]

mix_prob_prior <- c(0.5,0.5)
############## Initialize #####################
#spike_ind: binary indicator variable
error_sig2 <- 0.3^2

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

# x <- cbind(D_train,x_train)
# cv_model <- cv.glmnet(x, y, alpha = 1)
# best_lambda <- cv_model$lambda.min
# best_lambda
# best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
# coeff_lasso <- as.vector(coef(best_model))
# alpha <- coeff_lasso[2:(para+1)]
# df <- kmeans(alpha,3)
# alpha_group <- as.vector(df$centers)

#a1 <- rnorm(round(9*para/10),0,0.001)
a1 <- rep(0,round(9*para/10))
a2 <- rnorm(round(para/20),-0.5,0.01)
a3 <- rnorm(para-length(c(a1,a2)),0.5,0.01)

alpha <- sample(c(a1,a2,a3),para)

alpha_group <- c(0,-0.5,0.5)
mu <- 1 #prior mean of alpha
sigma2 <- 1 #prior variance of alpha

#
# reg <- lm(y~D_train)
# out <- reg$coefficients
# df <- kmeans(out,2)
#
# alpha <- out[-1]
# alpha_group <- as.vector(df$centers)
# mu <- alpha_group #prior mean of alpha
# sigma2 <- rep(0.1,ngroup) #prior variance of alpha

beta <- 1

prior_etax_var <- 0.3*diag(cov)
B=500
#=======================================

iteration <- 1500
burn <- 700

count <- 0

mix_prob_sim  <- array(0, c(ngroup,1, (iteration-burn)))
alpha_gr_sim <- array(0, c((ngroup+1), 1, (iteration-burn)))
alpha_sim <- array(0, c(para, 1, (iteration-burn)))
inv_img_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int1_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int2_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int3_sim <- array(0, c(dim, dim, (iteration-burn)))
inv_int4_sim <- array(0, c(dim, dim, (iteration-burn)))
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

q <- 0.5# This is important
err_alpha <- 1
err_beta <- 1
###################################################

for(ind in 1:iteration){
  alpha_img <- rep(0,dim^2)
  alpha_int1 <- rep(0,dim^2)
  alpha_int2 <- rep(0,dim^2)
  alpha_int3 <- rep(0,dim^2)
  alpha_int4 <- rep(0,dim^2)
  beta <- rnorm(1,(sum(y- (x_train %*% eta_x) - (D_train %*%alpha) )*prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)),(error_sig2 * prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)) )
  # beta <- rnorm(1,(sum(y- (x_train * eta_x) - (D_train %*%alpha) )*prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)),(error_sig2 * prior_beta_var)/(error_sig2+(ntrain *prior_beta_var)) )
  mm1 <- solve((t(x_train)%*% ((1/error_sig2)*diag(ntrain)) %*% x_train ) + solve(prior_etax_var)) %*%
    t(x_train)%*% ((1/error_sig2)*diag(ntrain))%*%(y - beta - (D_train %*%alpha) )
  mm2 <- solve((t(x_train)%*% ((1/error_sig2)*diag(ntrain)) %*% x_train ) + solve(prior_etax_var))
  eta_x <- mvrnorm(1,mm1,mm2)
  # for(h in 1:ngroup){
  #   sum1 <- NULL
  #
  #   for(k in 1:para){
  #
  #     # tmp1 <- numeric()
  #     # for(i in 1:ntrain){
  #     #   tmp1[i] <- sum(alpha*D_train[i,]) - (alpha[k]*D_train[i,k])
  #     # }
  #     tmp1 <- D_train%*%alpha - alpha[k]*D_train[,k]
  #
  #     sum1[k] <- mix_prob[h]*exp( (-1/(2*error_sig2)) *
  #                                   (((alpha_group[h]^2)*sum(D_train[,k]^2))- (2*alpha_group[h]*sum(D_train[,k]*(y- beta-(x_train %*% eta_x) -tmp1))) ) )
  #   }
  #   gr[,h] <- sum1
  # }
  
  for(k in 1:para){
    gr <- numeric()
    # c1 <- matrix(rep(D_train%*%alpha,para),nrow=ntrain)-hadamard.prod(matrix(rep(alpha,ntrain),nrow=ntrain,byrow = T),D_train)
    c2 <- sum(D_train[,k]^2)
    c3 <- sum(D_train[,k]*(y- beta-(x_train %*% eta_x)-((D_train%*%alpha)-(D_train[,k]*alpha[k]))))
    #c3 <- alpha_group[h+1]*apply(hadamard.prod(D_train,matrix(rep(y- beta-(x_train * eta_x),para),nrow=ntrain) -c1),2,sum  )
    gr <- (-1/(2*error_sig2)) *(((alpha_group[2:(ngroup+1)]^2)*c2)-(2*alpha_group[2:(ngroup+1)]*c3))
    gr <- c(0,gr)
    if (max(gr) < -750 | max(gr) > 700) gr <- gr + (700 - max(gr))  
    #gr1 <- t(as.matrix(apply(gr, 1, FUN=function(x) normalize(x) )))
    
    mix1 <- c((1-q),  q*mix_prob)
    
    prob <- hadamard.prod(exp(gr),mix1)
    
    class_alpha <- which(prob==max(prob))[1]
    
    alpha[k] <- alpha_group[class_alpha]
    #print(k)
  }
  
  
  q <- rbeta(1,(para-length(which(alpha==0))) + 1,length(which(alpha==0)) + 1 )
  
  # for(k in 1:para){
  #    class <- which(gr[k,]==max(gr[k,]))
  #    alpha[k] <- alpha_group[class]
  #  }
  
  temp1 <- numeric()
  for(h in 1:(ngroup)){
    temp1[h] <- length(which(alpha==alpha_group[h+1]))
  }
  
  mix_prob <- as.vector(rdirichlet(1,(temp1+mix_prob_prior)))
  
  for(h in 1:(ngroup)){
    if(temp1[h]>0){
      for(i in 1:ntrain){
        index <-  which(alpha==alpha_group[h+1])
        index_compli <-  which(alpha!=alpha_group[h+1])
        par1[i,h] <- sum(D_train[i,index])
        par1_compli[i,h] <- sum(alpha[index_compli]*D_train[i,index_compli])
      }
      
      a1 <- sum((par1[,h])^2)
      a2 <- sum(par1[,h]*(y-rep(beta,ntrain)- (x_train %*% eta_x)-par1_compli[,h] ))
      #a2 <- sum(par1[,h]*(y-rep(beta,ntrain)- (x_train * eta_x)-par1_compli[,h] ))
      mm <- ((a2/error_sig2) + (mu/sigma2))/((a1/error_sig2) + (1/sigma2))
      ss <- 1/((a1/error_sig2) + (1/sigma2))
      alpha_group[h+1] <- rnorm(1,mm,ss)
    }
    else{
      alpha_group[h+1] <- rnorm(1,mu, sigma2)
      
    }
  }
  mu <-  rnorm(1,solve(1+(length(unique(class_alpha))/sigma2))%*%(1+ sum(alpha_group[unique(class_alpha)]/sigma2)),solve(1+(length(unique(class_alpha))/sigma2)))
  sigma2 <- rinvgamma(1,shape=(length(unique(class_alpha))/2)+prior_sig_alpha,scale= (sum((alpha_group[unique(class_alpha)] - mu)^2)/2)+prior_sig_beta  )
  
  error_sig2 <- rinvgamma(1,shape=(ntrain/2)+err_alpha,scale= (sum((y-rep(beta,ntrain)- (x_train %*% eta_x) - (D_train%*%alpha))^2)/2 )+ err_beta)
  #error_sig2 <- rinvgamma(1,shape=(ntrain/2)+err_alpha,scale= (sum((y-rep(beta,ntrain)- (x_train * eta_x) - (D_train%*%alpha))^2)/2 )+ err_beta)
  alpha_img[ind_c1] <- alpha[1:wavelet]
  alpha_int1[ind_c1] <- alpha[(wavelet+1):(2*wavelet)]
  alpha_int2[ind_c1] <- alpha[((2*wavelet)+1):(3*wavelet)]
  alpha_int3[ind_c1] <- alpha[((3*wavelet)+1):(4*wavelet)]
  alpha_int4[ind_c1] <- alpha[((4*wavelet)+1):para]
  
  inv_img <- reconstr.2D(alpha_img,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int1 <- reconstr.2D(alpha_int1,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int2 <- reconstr.2D(alpha_int2,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int3 <- reconstr.2D(alpha_int3,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
  inv_int4 <- reconstr.2D(alpha_int4,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)#matrix.heatmap(inv_img)
  #matrix.heatmap(inv_img)
  
  if (ind > burn)
  {
    count <- count + 1
 
    mix_prob_sim[,,count]  <- mix_prob
    beta_sim[,,count] <- beta
    alpha_sim[,,count] <- alpha
    inv_img_sim [,,count] <- inv_img
    inv_int1_sim [,,count] <- inv_int1
    inv_int2_sim [,,count] <- inv_int2
    inv_int3_sim [,,count] <- inv_int3
    inv_int4_sim [,,count] <- inv_int4
    
    alpha_gr_sim[,,count] <- alpha_group
    eta_x_sim[,,count] <- eta_x
  }
  print(ind)
  print(table(alpha))
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
beta_est <- apply(beta_sim, 1, mean)
eta_x_est <- apply(eta_x_sim, 1, mean)

# Prediction error
y_pred <- numeric()
for(i in (ntrain+1):nind){
  y_pred[i-ntrain] <- beta_est + sum(D_full[i,]*alpha_est)+ sum(x_full[i,]*eta_x_est)
}
y_test <- y_full[(ntrain+1):nind]

pred_err_m1 <- sqrt(sum((y_test - y_pred)^2)) /sqrt(sum((y_test - mean(y_test))^2))
pred_err_m1

# Coverage
forecast <- forecast_cred(pred_sim,y_test)
cred <- forecast$cred_pred
covergae_prop <- length(which(cred[,3]<cred[,1] & cred[,1]<cred[,7]))/length(cred[,3])
width <- cred[,7]-cred[,3]
