#Generate square shaped functional regression coefficient
family = "DaubLeAsymm"
fil.num = 8
level = 3
dims=c(64,64)
p = prod(dims)

# generate coefficient images
u = seq(0,1,length.out=(dims[1]+1))[-1]
v = seq(0,1,length.out=(dims[2]+1))[-1]
beta1.0 = matrix(rep(0,prod(dims)),dims[1],dims[2])
for(a in 1:dims[1]){
  for(b in 1:dims[2]){
    if(u[a]>=(4*pi/80) & u[a]<=(14*pi/80) & v[b]>=(4*pi/80) & v[b]<=(14*pi/80)){
      beta1.0[a,b] = 3/5
    }
    if(u[a]>=(14*pi/80) & u[a]<=(24*pi/80) & v[b]>=(14*pi/80) & v[b]<=(24*pi/80)){
      beta1.0[a,b] = 1
    }
  }
}
matrix.heatmap(beta1.0)

wavelet1 = 64^2
family = "DaubLeAsymm"
fil.num = 8 #Can use 4/6/8/10
level = 3

beta_sq <- rep(wavelet1)
beta_sq <- decomp.2D(list(beta1.0),min.level=level,family=family,fil.num=fil.num)

tt <- reconstr.2D(beta_sq,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
matrix.heatmap(tt)

range(beta_sq)
length(beta_sq)
beta_sq[ which(beta_sq<(-0.5) )]=-3
beta_sq[ which(beta_sq>0.5)]= 3
beta_sq[ which(beta_sq !=(3) & beta_sq!=(-3))]=0

##Generate round shaped functional regression coefficient

family = "DaubLeAsymm"
fil.num = 8
level = 3

dims=c(64,64)
p = prod(dims)

# generate coefficient images
u = seq(0,1,length.out=(dims[1]+1))[-1]
v = seq(0,1,length.out=(dims[2]+1))[-1]

beta1.0 = matrix(rep(0,prod(dims)),dims[1],dims[2])
for(a in 1:dims[1]){
  for(b in 1:dims[2]){
    if((u[a]-9*pi/80)^2+(v[b]-9*pi/80)^2<=(3*pi/40)^2){
      beta1.0[a,b] = 3/5
    }
    if((u[a]-19*pi/80)^2+(v[b]-17*pi/80)^2<=(3*pi/40)^2){
      beta1.0[a,b] = 1
    }
  }
}

matrix.heatmap(beta1.0)

wavelet1 = 64^2
family = "DaubLeAsymm"
fil.num = 8 #Can use 4/6/8/10
level = 3

beta_sq <- rep(wavelet1)
beta_sq <- decomp.2D(list(beta1.0),min.level=level,family=family,fil.num=fil.num)

tt <- reconstr.2D(beta_sq,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
matrix.heatmap(tt)

range(beta_sq)
length(beta_sq)
beta_sq[ which(beta_sq<(-0.5) )]=-3
beta_sq[ which(beta_sq>0.5)]= 3
beta_sq[ which(beta_sq !=(3) & beta_sq!=(-3))]=0

##Generate triangle shaped functional regression coefficient
family = "DaubLeAsymm"
fil.num = 8
level = 3

dims=c(64,64)
p = prod(dims)

# generate coefficient images
u = seq(0,1,length.out=(dims[1]+1))[-1]
v = seq(0,1,length.out=(dims[2]+1))[-1]
# group 1
beta1.0 = matrix(rep(0,prod(dims)),dims[1],dims[2])
for(a in 1:dims[1]){
  for(b in 1:dims[2]){
    if(u[a]>=1/10 & u[a]<=5/10 & (u[a]/2+v[b])>=7/20 & (v[b]-u[a]/2)<=1/4){
      beta1.0[a,b] = 3/5
    }
    if(u[a]>=5/10 & u[a]<=9/10 & (u[a]/2+v[b])>=19/20 & (v[b]-u[a]/2)<=9/20){
      beta1.0[a,b] = 1
    }
  }
}

matrix.heatmap(beta1.0)

wavelet1 = 64^2
family = "DaubLeAsymm"
fil.num = 8 #Can use 4/6/8/10
level = 3

beta_sq <- rep(wavelet1)
beta_sq <- decomp.2D(list(beta1.0),min.level=level,family=family,fil.num=fil.num)

tt <- reconstr.2D(beta_sq,dims,family="DaubLeAsymm",fil.num=fil.num,min.level=level)
matrix.heatmap(tt)

beta_sq[ which(beta_sq<(-0.5) )]=-3
beta_sq[ which(beta_sq>0.7)]= 3
beta_sq[ which(beta_sq !=(3) & beta_sq!=(-3))]=0
