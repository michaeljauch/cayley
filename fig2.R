source('cayley_utilities.R')
library(MASS)
library(rstiefel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(5)

n <- 100 # Number of observations 
p <- 50 # Number of variables
k <- 3 # Dimension of low rank component of the covariance matrix
m <- p*k - k*(k+1)/2 # Dimension of the Stiefel manifold 

# Special matrices, computed using functions from cayley_utilities.R
Dk <- D_star(k) 
K <- as.matrix(K_rm(p-k, k))


# Set true values of parameters and generate a data matrix Y. 
Lambda <- diag(c(5,3,1.5))
Lambdainv <- diag(1/c(5,3,1.5))
sig2 <- 1
Q <- rustiefel(p,k)
Sigma <- sig2*(Q %*% Lambda %*% t(Q) + diag(p))
Y <- mvrnorm(n=n, mu=rep(0,p), Sigma = Sigma)

# Compute parameters of the matrix Bingham posterior. 
Abmf<- t(Y) %*% Y
Bbmf <- solve(Lambdainv + diag(k))/(2*sig2)
Cbmf <- matrix(rep(0, p*k), nrow=p)

# Calculate the posterior mode V. 
V <- eigen(Abmf)$vectors 


# Set the length/number of chains and decide how much burn-in to allow
num_chains <- 1
num_sim <- 10000
warmup <- 2000

# Simulate with Stan and extract simulated values
dat <- list('k'=k, 'p'=p, 'Dk'=Dk, 'K'=K, 'Abmf'=Abmf, 'Bbmf'=Bbmf, 'Cbmf'=Cbmf)
model_path <- 'stiefel_bmf.stan'
fit <- stan(file = model_path, data = dat, iter = num_sim + warmup, 
            chains = num_chains, refresh=1, warmup=warmup, seed=5)
draws <- extract(fit)
Q <- draws$Q

# Compute principal angles
angle_mat <- matrix(rep(0, num_sim*num_chains*k), nrow=num_sim*num_chains)
for(i in 1:(num_sim*num_chains)){
  for(j in 1:k){
    angle_mat[i,j] <- acos(abs(Q[i,,j] %*% V[,j])/
                             sqrt(Q[i,,j] %*% Q[i,,j])/
                             sqrt(V[,j] %*% V[,j]))
  }
}
angle_mat <- as.data.frame(angle_mat)
names(angle_mat) <- c('theta1', 'theta2', 'theta3')


# Now simulate from the same distribution using the rstiefel package
Qhoff <- array(rep(0,(num_sim + warmup)*p*k), dim=c(num_sim + warmup, p, k))
Qhoff[1,,] <- rustiefel(p,k)
for(t in 2:(num_sim + warmup)){
  print(t)
  Qhoff[t,,] <- rbmf.matrix.gibbs(A=Abmf, B=Bbmf, C=Cbmf, X=Qhoff[t-1,,])
}
Qhoff <- Qhoff[(warmup+1):(num_sim + warmup),,]

# Compute principal angles
angle_mat_hoff <- matrix(rep(0, num_sim*k), nrow=num_sim)
for(i in 1:num_sim){
  for(j in 1:k){
    angle_mat_hoff[i,j] <- acos(abs(Qhoff[i,,j] %*% V[,j])/
                                  sqrt(Qhoff[i,,j] %*% Qhoff[i,,j])/
                                  sqrt(V[,j] %*% V[,j]))
  }
}
angle_mat_hoff <- as.data.frame(angle_mat_hoff)
names(angle_mat_hoff) <- c('theta1', 'theta2', 'theta3')

# Begin creating Figure
pdf(file='Binghamsim.pdf', width=9, height=7, family='Times')
arrangement <- rbind(c(1,1,1,2,2,2), c(3,3,4,4,5,5))
layout(arrangement)

# Trace plot
par(cex=1.3, mai=c(1,1,.5,.5), mgp=c(1.75, .75, 0))
plot(x=angle_mat[1:2000,1], type='l', lwd=1.5, ylab=expression(theta[1]), 
     xlab='Index')
points(angle_mat_hoff[1:2000,1], type='l', ylim=c(0,.5), col="grey75", lwd=2)

# ACF plot
par(mai=c(1,1,.5,.5), mgp=c(1.75, .75, 0))
acf_obj <- acf(angle_mat[,1], plot=FALSE)
acf_obj_hoff <- acf(angle_mat_hoff[,1], plot=FALSE)
plot(acf_obj$lag, acf_obj$acf, ylim=c(-.1, 1), xlab='Lag', 
     ylab='Autocorrelation', pch=19)
points(acf_obj_hoff$lag, acf_obj_hoff$acf, ylim=c(-.1, 1), col="grey75", pch=19)

# Margin plots
par(mai=c(1,1,.25,.5), mgp=c(1.75, .75, 0))
mintheta1 <- min(angle_mat$theta1)
maxtheta1 <- max(angle_mat$theta1)
hist(~theta1,data=angle_mat, col=rgb(0.1,0.1,0.1,0.5), main='', 
     xlab=expression(theta[1]), ylab='Density', 
     freq=FALSE, breaks=seq(mintheta1*.8,maxtheta1*1.2, length.out=15), 
     xlim=c(mintheta1*.8,maxtheta1*1.2))
hist(~theta1, data=angle_mat_hoff, col=rgb(0.8,0.8,0.8,0.5), add=T, freq=FALSE, 
     ylab='', breaks=seq(mintheta1*.8,maxtheta1*1.2, length.out=15), 
     xlim=c(mintheta1*.8,maxtheta1*1.2))

par(mai=c(1,1,.25,.5), mgp=c(1.75, .75, 0))
mintheta2 <- min(angle_mat$theta2)
maxtheta2 <- max(angle_mat$theta2)
hist(~theta2, data=angle_mat, col=rgb(0.1,0.1,0.1,0.5), main='', 
     xlab=expression(theta[2]), ylab='', freq=FALSE, 
     breaks=seq(mintheta2*.8,maxtheta2*1.2, length.out=15), 
     xlim=c(mintheta2*.8,maxtheta2*1.2))
hist(~theta2, data=angle_mat_hoff, col=rgb(0.8,0.8,0.8,0.5), add=T, freq=FALSE, 
     ylab='',  breaks=seq(mintheta2*.8,maxtheta2*1.2, length.out=15), 
     xlim=c(mintheta2*.8,maxtheta2*1.2))

par(mai=c(1,1,.25,.5), mgp=c(1.75, .75, 0))
mintheta3 <- min(angle_mat$theta3)
maxtheta3 <- max(angle_mat$theta3)
hist(~theta3, data=angle_mat, col=rgb(0.1,0.1,0.1,0.5), main='', 
     xlab=expression(theta[3]), ylab='', freq=FALSE, 
     breaks=seq(mintheta3*.8,maxtheta3*1.2, length.out=15), 
     xlim=c(mintheta3*.8,maxtheta3*1.2))
hist(~theta3, data=angle_mat_hoff, col=rgb(0.8,0.8,0.8,0.5), add=T, freq=FALSE, 
     breaks=seq(mintheta3*.8,maxtheta3*1.2, length.out=15), 
     xlim=c(mintheta3*.8,maxtheta3*1.2), ylab='')
dev.off()


