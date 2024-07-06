library(oro.nifti)
library(neurobase)
library(plsgenomics)
library(OpenImageR)

n=71
load("image_CT.dat")
img_1 <- array(0,c(64,64,n))
img_2 <- array(0,c(64,64,n))
img_3 <- array(0,c(64,64,n))
img_4 <- array(0,c(64,64,n))
for(i in 1:n){
  img_1[,,i] <- image[[i]][1:64,1:64]
}
for(i in 1:n){
  img_2[,,i] <- image[[i]][1:64,65:128]
}
for(i in 1:n){
  img_3[,,i] <- image[[i]][65:128,1:64]
}
for(i in 1:n){
  img_4[,,i] <- image[[i]][65:128,65:128]
}
mean_1=apply(img_1,c(1,2),mean)
sd_1=apply(img_1,c(1,2),sd)
sd_1[which(sd_1==0 )]=1

mean_2=apply(img_2,c(1,2),mean)
sd_2=apply(img_2,c(1,2),sd)
sd_2[which(sd_2==0 )]=1

mean_3=apply(img_3,c(1,2),mean)
sd_3=apply(img_3,c(1,2),sd)
sd_3[which(sd_3==0 )]=1

mean_4=apply(img_4,c(1,2),mean)
sd_4=apply(img_4,c(1,2),sd)
sd_4[which(sd_4==0 )]=1

##Normalize data

for(i in 1:n){
  img_1[,,i] <- (img_1[,,i] -  mean_1)/sd_1
}
for(i in 1:n){
  img_2[,,i] <- (img_2[,,i] -  mean_2)/sd_2
}
for(i in 1:n){
  img_3[,,i] <- (img_3[,,i] -  mean_3)/sd_3
}
for(i in 1:n){
  img_4[,,i] <- (img_4[,,i] -  mean_4)/sd_4
}

wavelet = 64^2
nind=n
family = "DaubLeAsymm"
fil.num = 8 #Can use 4/6/8/10
level = 3

wv_1 <- matrix(0,nrow=n,ncol=wavelet)
for(i in 1:nind)
{
  wv_1[i,] <- decomp.2D(list(img_1[,,i]),min.level=level,family=family,fil.num=fil.num)
}

dm_wv1 <- dim(wv_1)[2]

tt <- abs(apply(wv_1,2,sd))
ind2 <- (which(tt<0.08))
ind_c1 <- c(1:dm_wv1)[-ind2]
length(ind_c1)#446
wv1 <- wv_1[,ind_c1]
dm_wv2 <- dim(wv1)[2]
dm_wv2

save(wv1,file="wv_part1.dat")
save(ind_c1,file="ind_c1_part1.dat")

#########################


wv_2 <- matrix(0,nrow=n,ncol=wavelet)
for(i in 1:nind)
{
  wv_2[i,] <- decomp.2D(list(img_2[,,i]),min.level=level,family=family,fil.num=fil.num)
}

dm_wv1 <- dim(wv_2)[2]
tt <- abs(apply(wv_2,2,sd))
ind2 <- (which(tt<0.07))
ind_c1 <- c(1:dm_wv1)[-ind2]

wv1 <- wv_2[,ind_c1]
dm_wv2 <- dim(wv1)[2]
dm_wv2

save(wv1,file="wv_part2.dat")
save(ind_c1,file="ind_c1_part2.dat")

################################

wv_3 <- matrix(0,nrow=n,ncol=wavelet)

for(i in 1:nind)
{
  wv_3[i,] <- decomp.2D(list(img_3[,,i]),min.level=level,family=family,fil.num=fil.num)
}

dm_wv1 <- dim(wv_3)[2]

tt <- abs(apply(wv_3,2,sd))
ind2 <- (which(tt<0.09))
ind_c1 <- c(1:dm_wv1)[-ind2]
wv1 <- wv_3[,ind_c1]
dm_wv2 <- dim(wv1)[2]
dm_wv2
save(wv1,file="wv_part3.dat")
save(ind_c1,file="ind_c1_part3.dat")

#######################


wv_4 <- matrix(0,nrow=n,ncol=wavelet)

for(i in 1:nind)
{
  wv_4[i,] <- decomp.2D(list(img_4[,,i]),min.level=level,family=family,fil.num=fil.num)
}
dm_wv1 <- dim(wv_4)[2]

tt <- abs(apply(wv_4,2,sd))
ind2 <- (which(tt<0.07))
ind_c1 <- c(1:dm_wv1)[-ind2]

wv1 <- wv_4[,ind_c1]
dm_wv2 <- dim(wv1)[2]
dm_wv2

save(wv1,file="wv_part4.dat")
save(ind_c1,file="ind_c1_part4.dat")
