A <- as.matrix(kronecker(rep(1, n), alpha))
B <- as.matrix(kronecker(rep(1, n), beta))
C <- as.matrix(kronecker(theta, gamma))

### estimation ----
#### FF
at <- Rt <- mt <- Ct <- rep(NA, t)
ft <- matrix(NA, ncol = t, nrow = ages*n)
Qt <- array(NA, dim=c(n*ages, n*ages, t))

at[1] <- nu + m0
Rt[1] <- C0 + sigma_w ##phiw

ft[,1] <- A + (B * at[1]) + C
#Qt[,,1] <- Rt[1] * B %*% t(B) + diag(n*ages)*sigma_e #phie
Qt[,,1] <- Rt[1] * tcrossprod(B, B) + diag(n*ages)*sigma_e #phie

#mt[1] <- at[1] + Rt[1] * ( t(B) %*% chol2inv(chol(Qt[,,1])) %*% (y_obs[,1] - ft[,1]) )
mt[1] <- at[1] + Rt[1] * ( crossprod(B, chol2inv(chol(Qt[,,1]))) %*% (y_obs[,1] - ft[,1]) )
#Ct[1] <- Rt[1] - Rt[1]^2 * ( t(B) %*% chol2inv(chol(Qt[,,1])) %*% B )
Ct[1] <- Rt[1] - Rt[1]^2 * ( crossprod(B, chol2inv(chol(Qt[,,1]))) %*% B )

t = 20

ff1 <- function(t){
  for(k in 2:t){
    at[k] <- nu + mt[k-1]
    Rt[k] <- Ct[k-1] + sigma_w ##phiw
    
    ft[,k] <- A + (B * at[k]) + C
    Qt[,,k] <- Rt[k] * tcrossprod(B) + diag(n*ages)*sigma_e #phie
    
    mt[k] <- at[k] + Rt[k] * ( crossprod(B, chol2inv(chol(Qt[,,k]))) %*% (y_obs[,k] - ft[,k]) )
    Ct[k] <- Rt[k] - Rt[k]^2 * ( crossprod(B, chol2inv(chol(Qt[,,k]))) %*% B )
  }
}

ff2 <- function(t){
  for(k in 2:t){
    at[k] <- nu + mt[k-1]
    Rt[k] <- Ct[k-1] + sigma_w ##phiw
    
    ft[,k] <- A + (B * at[k]) + C
    Qt[,,k] <- Rt[k] * B %*% t(B) + diag(n*ages)*sigma_e #phie
    
    mt[k] <- at[k] + Rt[k] * ( t(B) %*% chol2inv(chol(Qt[,,k])) %*% (y_obs[,k] - ft[,k]) )
    Ct[k] <- Rt[k] - Rt[k]^2 * ( t(B) %*% chol2inv(chol(Qt[,,k])) %*% B )
  }
}

ff3 <- function(t){
  for(k in 2:t){
    at[k] <- nu + mt[k-1]
    Rt[k] <- Ct[k-1] + sigma_w ##phiw
    
    ft[,k] <- A + (B * at[k]) + C
    Qt[,,k] <- Rt[k] * tcrossprod(B) + diag(n*ages)*sigma_e #phie
    
    mt[k] <- at[k] + Rt[k] * tcrossprod( crossprod(B, chol2inv(chol(Qt[,,k]))), t(y_obs[,k] - ft[,k]) )
    Ct[k] <- Rt[k] - Rt[k]^2 * tcrossprod( crossprod(B, chol2inv(chol(Qt[,,k]))), t(B) )
  }
}


library(microbenchmark)
microbenchmark(ff1(t),
               ff2(t),
               ff3(t),
               times = 100)

### FF1 mais rapido, seguido por FF3 e FF2 na media
### FF1, FF2 e FF3 na mediana


### TESTAR PACOTE FKF (fast kalman filter)
library(FKF)
