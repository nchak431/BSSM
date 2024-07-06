decomp.2D <- function(x,family="DaubLeAsymm",fil.num=4,min.level=2){
  # if((class(x)!="list")||(class(x[[1]])!="matrix")){
  #   cat("input should be a list of square matrices","\n")
  # }
  nsample = length(x)
  dims = dim(x[[1]])
  upper.level = min(log2(dims))
  if(upper.level<=min.level){
    cat("invalid min.level, should be less than min(log2(dim(x[[1]]))).","\n")
  }
  dim1d = prod(dims)
  wdt = matrix(NA,nsample,dim1d)

  for(i in 1:nsample){
    input = x[[i]]
    wd0 = imwd(input,filter.number=fil.num,family=family)

    for(j in (upper.level-1):min.level){
      rep.sub = c(eval(parse(text=paste0("wd0$w",j,"L4"))),eval(parse(text=paste0("wd0$w",j,"L1"))),
                  eval(parse(text=paste0("wd0$w",j,"L2"))),eval(parse(text=paste0("wd0$w",j,"L3"))))
      wdt[i,1:(4^(j+1))] = rep.sub
    }

  }

  return(wdt)
}

reconstr.2D <- function(eta,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level){
  if(length(eta)!=prod(dims)){cat("length of eta does not match the given dims.","/n")}
  upper.level = min(log2(dims))
  wr0 = imwd(matrix(rep(0,prod(dims)),dims[1],dims[2]),filter.number=fil.num,family=family,RetFather=F)
  
  if(min.level>1){
    sub.wd = imwd(matrix(eta[1:(4^min.level)],2^min.level,2^min.level),family=family,filter.number=fil.num)
    wr0$w0Lconstant = sub.wd$w0Lconstant
    for(j in 0:(min.level-1)){
      value = eval(parse(text=paste0("sub.wd$w",j,"L1")))
      eval(parse(text=paste0("wr0$w",j,"L1 = value")))
      value = eval(parse(text=paste0("sub.wd$w",j,"L2")))
      eval(parse(text=paste0("wr0$w",j,"L2 = value")))
      value = eval(parse(text=paste0("sub.wd$w",j,"L3")))
      eval(parse(text=paste0("wr0$w",j,"L3 = value")))
    }
    for(j in min.level:(upper.level-1)){
      size = 4^j
      eval(parse(text=paste0("wr0$w",j,"L1 = eta[(size+1):(2*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L2 = eta[(2*size+1):(3*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L3 = eta[(3*size+1):(4*size)]")))
    }
    wr = imwr(wr0)
    
  }else if(min.level==1){
    a = eta[1]; b = eta[2]; c = eta[3]; d = eta[4]
    wr0$w0Lconstant = wr0$w0L4 = (a+b+c+d)/2
    wr0$w0L1 = (a+b-c-d)/2
    wr0$w0L2 = (a+c-b-d)/2
    wr0$w0L3 = (a+d-b-c)/2
    for(j in 1:(upper.level-1)){
      size = 4^j
      eval(parse(text=paste0("wr0$w",j,"L1 = eta[(size+1):(2*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L2 = eta[(2*size+1):(3*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L3 = eta[(3*size+1):(4*size)]")))
    }
    wr = imwr(wr0)
    
  }else{
    for(j in 0:(upper.level-1)){
      size = 4^j
      eval(parse(text=paste0("wr0$w",j,"L1 = eta[(size+1):(2*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L2 = eta[(2*size+1):(3*size)]")))
      eval(parse(text=paste0("wr0$w",j,"L3 = eta[(3*size+1):(4*size)]")))
    }
    wr = imwr(wr0)
  }
  
  return(wr)
}
forecast_cred <- function(sim,true)
{
  sim_q11 <- apply(sim,c(1,2),quantile,probs=0.025)
  sim_q1 <- apply(sim,c(1,2),quantile,probs=0.25)
  sim_med <- apply(sim,c(1,2),quantile,probs=0.50)
  sim_q3 <- apply(sim,c(1,2),quantile,probs=0.75)
  sim_q33 <- apply(sim,c(1,2),quantile,probs=0.975)
  sim_mean <- apply(sim,c(1,2),mean)
  sim_l <- apply(sim,c(1,2),quantile,probs=0.005)
  sim_u <- apply(sim,c(1,2),quantile,probs=0.995)
  
  err <-  (sum((true-sim_med)^2))/(sum((true - mean(true))^2))
  err1 <- sqrt(sum((true-sim_med)^2))/sqrt(sum((true - mean(true))^2))
  
  cred <- cbind(true,sim_med,sim_q11,sim_q1,sim_med,sim_q3,sim_q33,sim_l,sim_u,sim_mean)
  colnames(cred)=c("true","est","2.5%","25%","50%","75%","97.5%","0.5","99.5","mean")
  
  return(list(med_est = sim_med,pred_err=err,sq_pred_err=err1,cred_pred=cred ))
}

cred_voxel <- function(sim)
{
  sim_q11 <- apply(sim,c(1,2),quantile,probs=0.025)
  sim_q1 <- apply(sim,c(1,2),quantile,probs=0.25)
  sim_med <- apply(sim,c(1,2),quantile,probs=0.50)
  sim_q3 <- apply(sim,c(1,2),quantile,probs=0.75)
  sim_q33 <- apply(sim,c(1,2),quantile,probs=0.975)
  sim_mean <- apply(sim,c(1,2),mean)
  
  q11 <- vec(sim_q11);q1 <- vec(sim_q1);med <- vec(sim_med);q3 <- vec(sim_q3);q33 <- vec(sim_q33);mean <- vec(sim_mean)
  
  cred <- cbind(q11,q1,med,q3,q33,mean)
  colnames(cred)=c("2.5%","25%","50%","75%","97.5%","mean")
  
  return(list(cred=cred ))
}
alpha_draw <- function(alpha,beta,eta,sigma,alpha_group,q,mix_prob){
  for(k in 1:para){
    gr <- numeric()
    # c1 <- matrix(rep(D_train%*%alpha,para),nrow=ntrain)-hadamard.prod(matrix(rep(alpha,ntrain),nrow=ntrain,byrow = T),D_train)
    c2 <- sum(D_train[,k]^2)
    c3 <- sum(D_train[,k]*(y- beta-(x_train %*% eta)-((D_train%*%alpha)-(D_train[,k]*alpha[k]))))
    #c3 <- alpha_group[h+1]*apply(hadamard.prod(D_train,matrix(rep(y- beta-(x_train * eta_x),para),nrow=ntrain) -c1),2,sum  )
    gr <- (-1/(2*sigma)) *(((alpha_group[2:(ngroup+1)]^2)*c2)-(2*alpha_group[2:(ngroup+1)]*c3))
    gr <- c(0,gr)
    if (max(gr) < -750 | max(gr) > 700) gr <- gr + (700 - max(gr))  
    #gr1 <- t(as.matrix(apply(gr, 1, FUN=function(x) normalize(x) )))
    
    mix1 <- c((1-q),  q*mix_prob)
    
    prob <- hadamard.prod(exp(gr),mix1)
    
    class_alpha <- which(prob==max(prob))[1]
    
    alpha[k] <- alpha_group[class_alpha]
    #print(k)
  }
  alpha
}

clus_draw <- function(alpha,alpha_group,beta,eta,sigma,mu,sig){
for(h in 1:(ngroup)){
  if(temp1[h]>0){
    for(i in 1:ntrain){
      index <-  which(alpha==alpha_group[h+1])
      index_compli <-  which(alpha!=alpha_group[h+1])
      par1[i,h] <- sum(D_train[i,index])
      par1_compli[i,h] <- sum(alpha[index_compli]*D_train[i,index_compli])
    }
    
    a1 <- sum((par1[,h])^2)
    a2 <- sum(par1[,h]*(y-rep(beta,ntrain)- (x_train %*% eta)-par1_compli[,h] ))
    #a2 <- sum(par1[,h]*(y-rep(beta,ntrain)- (x_train * eta_x)-par1_compli[,h] ))
    mm <- ((a2/sigma) + (mu/sig))/((a1/sigma) + (1/sig))
    ss <- 1/((a1/sigma) + (1/sig))
    alpha_group[h+1] <- rnorm(1,mm,ss)
  }
  else{
    alpha_group[h+1] <- rnorm(1,mu, sigma2)
    
  }
}
  alpha_group
}
