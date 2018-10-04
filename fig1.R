source('cayley_utilities.R')
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(5)

k <- 3 # Number of columns of our orthogonal matrix

# Setting the length/number of chains and deciding how much burn-in to allow
num_chains <- 4
num_samples <- 2500
warmup <- 2000

# This function evaluates the density associated with the marginal distribution 
# of an entry of a uniform random orthogonal matrix. The density is given in 
# Eaton (1989), referenced in the paper. 
pdf_q11 <- function(q11, p){
  value <- 
    (2*pi)^(-1/2)*wishart_constant(p-1,1)/
    wishart_constant(p,1)*(1-q11^2)^((p-3)/2)
  idx <- (q11 <= 0 || q11 >= 1)*1
  value[idx] <- 0
  return(value)
}

# Calls a function from the file cayley_utilities.R which constructs the special 
# D matrix discussed in the appendix of the paper.
Dk <- D_star(k) 

pdf(file='~/Desktop/Unifsim.pdf', width=7, height=5, family='Times')
par(mfrow=c(2,2), cex=1, mgp=c(1.75, .75, 0), mai=c())

for(p in c(5, 50)){
  
  m <- p*k - k*(k+1)/2 # The dimension of the Stiefel manifold
  K <- as.matrix(K_rm(p-k, k)) # Constructs the commutation matrix
  dat <- list('k'=k, 'p'=p, 'Dk'=Dk, 'K'=K)
  model_path <- 'stiefel_unif.stan'
  fit <- stan(file = model_path, data = dat, iter = num_samples + warmup, 
              chains = 4, refresh=1, warmup=warmup, seed=5)
  draws <- extract(fit)
  
  # q11 plot + exact marginal pdf of the top left entry of Q
  x <- seq(-1,1,.01)
  density_q <- pdf_q11(x, p)
  Q <- draws$Q
  hist(c(Q[,1,1], -1, 1), xlim=c(-1,1), freq=FALSE, 
       main=paste('p =', p, ', k=3'), xlab=expression('Top left entry of'~Q), 
       breaks=25, ylim=c(0, max(density_q)*1.2))
  points(x, density_q, type='l', lwd=2)

  # phi plot, discarding simulated values far in the tails for the sake of 
  # having a sane plot. 
  phi <- draws$phi
  idx <- (sqrt(p/2)*phi[,1] >=-5 & sqrt(p/2)*phi[,1] <= 5)
  phi_without_tails <- phi[idx,1]
  x <- seq(-5, 5, .01)
  density_phi <- dnorm(x, mean=0, sd=1)
  hist(sqrt(p/2)*phi_without_tails, 
     xlab=expression('Rescaled first element of' ~varphi), breaks=25,
     freq=FALSE, main=paste('p =', p, ', k=3'),
     ylim=c(0, max(density_phi)*1.2), xlim=c(-5,5))
  points(x, density_phi, type='l', lwd=2)
}
dev.off()

