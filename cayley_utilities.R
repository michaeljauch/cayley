library(Matrix)

# This is w(r, p) from pg. 175 of Eaton (1983).

wishart_constant <- function(r,p){
  crp <- 2^(r*p/2)*pi^(p*(p-1)/4)
  for(j in 1:p){
    crp <- crp*gamma((r-j+1)/2)
  }
  return(1/crp)
}

# The following functions are used in the construction of the D matrix described
# in the appendix of the paper. The names of the functions come from Neudecker
# (1983). The notation does not align with that in the appendix. 

E <- function(i,j,k){
  E <- matrix(rep(0,k^2), nrow=k)
  E[i,j] <- 1
  return(E)
}

nu <- function(i,j,k){
  nu_v <- matrix(rep(0,k*(k-1)/2), ncol=1)
  nu_v[(j-1)*k + i -j*(j+1)/2,1] <- 1 
  return(nu_v)
}

D_star <- function(k){
  tmp <- matrix(rep(0,k^2*k*(k-1)/2), nrow=k^2)
  for(i in 1:k){
    for(j in 1:(i-1)){
      Ediff <- E(i,j,k) - E(j,i,k)
      vecEdiff <- as.matrix(c(Ediff), ncol=1)
      tmp <- tmp + vecEdiff %*% t(nu(i,j,k))
    }
  }
  return(tmp)
}

# The following functions are used in the construction of the commutation 
# matrix. The notation does not align with that in the appendix; however, this 
# is equivalent to the explicit construction given on pg. 37 of the book 'Linear 
# Structures' [Magnus (1988)]. 

e_rj <- function(r,j){
  tmp <- Matrix(rep(0,r), ncol=1, sparse=TRUE)
  tmp[j] <- 1
  return(tmp)
}

# Commutation matrix 
K_rm <- function(r,m){
  K_tmp <- Matrix(rep(0, r^2*m^2), nrow=r*m, sparse=TRUE)
  for(i in 1:r){
    for(j in 1:m){
      K_tmp <- K_tmp + (e_rj(r,i) %*% t(e_rj(m,j))) %x% 
        (e_rj(m,j) %*% t(e_rj(r,i)))
    }
  }
  return(K_tmp)
}