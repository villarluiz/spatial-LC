

ff_sp <- function(y_obs, ages, t, n, m0, C0, sigma_e, sigma_w, nu, alpha, beta, gamma, theta){
  ### aux ----
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
  Qt[,,1] <- Rt[1] * tcrossprod(B) + diag(n*ages)*sigma_e #phie
  
  #mt[1] <- at[1] + Rt[1] * ( t(B) %*% chol2inv(chol(Qt[,,1])) %*% (y_obs[,1] - ft[,1]) )
  aux <- t(B) %*% chol2inv(chol(Qt[,,1]))
  mt[1] <- at[1] + Rt[1] * ( aux %*% (y_obs[,1] - ft[,1]) )
  #Ct[1] <- Rt[1] - Rt[1]^2 * ( t(B) %*% chol2inv(chol(Qt[,,1])) %*% B )
  Ct[1] <- Rt[1] - Rt[1]^2 * ( aux %*% B )
  
  ##### FF ta demorando o tempo todo das iteracoes
  start <- proc.time()[[3]]
  for(k in 2:t){
    at[k] <- nu + mt[k-1]
    Rt[k] <- Ct[k-1] + sigma_w ##phiw
    
    ft[,k] <- A + (B * at[k]) + C
    #Qt[,,k] <- Rt[k] * B %*% t(B) + diag(n*ages)*sigma_e #phie
    Qt[,,k] <- Rt[k] * tcrossprod(B) + diag(n*ages)*sigma_e #phie
    
    #mt[k] <- at[k] + Rt[k] * ( t(B) %*% chol2inv(chol(Qt[,,k])) %*% (y_obs[,k] - ft[,k]) )
    aux <- t(B) %*% chol2inv(chol(Qt[,,k]))
    mt[k] <- at[k] + Rt[k] * ( aux %*% (y_obs[,k] - ft[,k]) )
    #Ct[k] <- Rt[k] - Rt[k]^2 * ( t(B) %*% chol2inv(chol(Qt[,,k])) %*% B )
    Ct[k] <- Rt[k] - Rt[k]^2 * ( aux %*% B )
  }
  end <- proc.time()[[3]]
  end - start
  return(list(mt = mt,
              Ct = Ct,
              Rt = Rt,
              at = at,
              ft = ft,
              Qt = Qt))
}

bs_sp <- function(mt, Ct, at, Rt, m0, C0){
  t = length(mt)
  kappa.bs <- rep(NA, t)
  kappa.bs[t] <- rnorm(1, mt[t], sqrt(Ct[t])) 
  
  #### BS
  st <- St <- rep(NA, t)
  st[t] <- mt[t]
  St[t] <- Ct[t]
  
  for(k in (t-1):1){
    Bt = Ct[k] * 1/Rt[k+1]
    st[k] = mt[k] + Bt %*% (kappa.bs[k+1] - at[k+1])
    St[k] = Ct[k] - Bt * Rt[k+1] * t(Bt)
    
    kappa.bs[k] <- rnorm(1, st[k], sqrt(St[k]))
  }
  
  Bt = C0 * 1/Rt[1]
  s0 <- m0 + Bt %*% (kappa.bs[1] - at[1])
  S0 <- C0 - Bt * Rt[1] * t(Bt)
  
  kappa0.bs <- rnorm(1, s0, sqrt(S0))
  
  return(list(kappa.bs = kappa.bs,
              kappa0.bs = kappa0.bs,
              st = st,
              St = St))
}

ffbs_sp <- function(y_obs, ages, t, n, m0, C0, sigma_e, sigma_w, nu, alpha, beta, gamma, theta){
  aux.ff <- ff_sp(y_obs, ages, t, n, m0, C0, sigma_e, sigma_w, nu, alpha, beta, gamma, theta)
  aux.bs <- bs_sp(aux.ff$mt, aux.ff$Ct, aux.ff$at, aux.ff$Rt, m0, C0)
  return(aux.bs)
}

#### SIMULACAO ----
library(sf)
library(spdep)
library(tidyverse)

### theta ----
rj_micro <- st_read("RJ_Microrregioes_2022/RJ_Microrregioes_2022.shp")
shape_aisp_nb_val <- st_make_valid(rj_micro[,c(2,5)])
nb <- poly2nb(shape_aisp_nb_val, queen=TRUE)

# matriz de vizinhaça W e M
W <- nb2mat(nb, style="B", zero.policy = TRUE)
M <- diag(rowSums(W))
# labda.dom <- eigen(W) #### contem 1, por ex.
# c(1/min(labda.dom$values), 1/max(labda.dom$values))
# set.seed(1)
#lambda <- runif(1, 1/min(labda.dom$values), 1/max(labda.dom$values))
#lambda = 0.15
#Q = M - lambda*W

ages <- 13  ## idades
t <- 40  ## tempo
n <- 15   ## regioes

Wi <- W[1:n, 1:n]
Mi <- diag(rowSums(Wi))
#labda.dom.i <- eigen(Wi) 
#c(1/min(labda.dom.i$values), 1/max(labda.dom.i$values)) ## ok

Qi = Mi - Wi
Qii = Qi[1:(n-1), 1:(n-1)]

sigma_t = 6
set.seed(14)
theta <- MASS::mvrnorm(n = 1, mu = rep(0, n-1), Sigma = chol2inv(chol(Qii))*sigma_t) 
theta_n <- rnorm(1, sum( (Wi[n, -n]/Mi[n,n])*theta[-n] ), sqrt(sigma_t/Mi[n,n]))
theta <- c(theta, theta_n)
### alpha, beta and kappa ----
## acho que tem constraint demais...
set.seed(14)
alpha <- readRDS("sim_ab.RDS")$alp
# #alpha <- alpha/sum(alpha)
#beta <- runif(ages, 0, 0.1)
beta <- readRDS("sim_ab.RDS")$bet
gamma <- runif(ages, 0, 0.1)
#beta <- beta/sum(beta)



kappa0 <- 1
kappa <- rep(NA, t)
nu <- -0.95

set.seed(14)
kappa[1] <- nu + kappa0 + rnorm(1, 0, sqrt(0.35)) 
for(k in 2:t){
  kappa[k] <- nu + kappa[k-1] + rnorm(1, 0, sqrt(0.35))
}
inc <- (kappa[1] - kappa[t]) / (t-1)

#alpha <- alpha + beta*mean(kappa) + gamma*mean(theta)
theta <- (theta - mean(theta))*sum(gamma)
#inclination <- (kappa[1] - kappa[t])/(t-1)
kappa <- (kappa - mean(kappa))/inc
#kappa <- aux.k[2:(t+1)]; kappa0 <- aux.k[1]
#beta <- beta*inc
gamma <- gamma/sum(gamma)


m0 <- 0; C0 <- 100
#sigma_e = runif(ages, min = 0, max = 1)
sigma_e = 0.2
sigma_w = 4

### aux ----
A <- kronecker(rep(1, n), alpha)
B <- kronecker(rep(1, n), beta)
C <- kronecker(theta, gamma)
set.seed(14)
y_obs <- matrix(NA, nrow = ages*n, ncol = t)
for(k in 1:t){
  y_obs[,k] <- A + (B * kappa[k]) + C + rnorm(n*ages, 0, sd = sqrt(sigma_e))
}

##

### estimation (1-step) #kappa ----
start <- proc.time()[[3]]
fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
               sigma_e = sigma_e, sigma_w = sigma_w, nu = nu, alpha = alpha,
               beta = beta, gamma = gamma, theta = theta)
end <- proc.time()[[3]]
end - start ## seconds

plot(1:t, kappa)
lines(1:t, fit$kappa.bs)
aux <- c(fit$kappa0.bs, fit$kappa.bs)
aux <- aux - mean(aux)
(aux[t+1] - aux[1])/t

### estimation (n-steps) #age-time-spatial parameters (no filtering, gamma + theta constraint)----
it = 20000 #iterations
kappa.chain <- matrix(NA, nrow = it, ncol = t)

alpha.chain <- matrix(NA, nrow = it, ncol = ages)
beta.chain <- matrix(NA, nrow = it, ncol = ages)
gamma.chain <- matrix(NA, nrow = it, ncol = ages)
sigma_e.chain <- matrix(NA, nrow = it, ncol = ages)

theta.chain <- matrix(NA, nrow = it, ncol = n)

sigma_w.chain <- rep(NA, it)
sigma_t.chain <- rep(NA, it)
#lambda.chain <- rep(NA, it)
nu.chain <- rep(NA, it)

## starters
alpha.chain[1, ] <- runif(ages)
beta.chain[1, ] <- runif(ages)
gamma.chain[1, ] <- runif(ages)

theta.chain[1, ] <- runif(n)    ### theta = 1 1 1 causou instabilidade

sigma_w.chain[1] <- 1
sigma_e.chain[1, ] <- rep(1, ages)
sigma_t.chain[1] <- 1

nu.chain[1] <- runif(1, -1, 1)

## gibbs
pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = it, clear = FALSE, width = 60)
pb$tick()
prop = 0

set.seed(16)
for(i in 2:it){
  
  pb$tick()
  ### kappa
  # fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
  #                sigma_e = sigma_e.chain[i-1],
  #                sigma_w = sigma_w.chain[i-1],
  #                nu = nu.chain[i-1],
  #                alpha = alpha.chain[i-1, ],
  #                beta = beta.chain[i-1, ],
  #                gamma = gamma.chain[i-1, ],
  #                theta = theta.chain[i-1, ])
  # 
  # kappa.chain[i,] <- fit$kappa.bs
  ### kappa constraint
  #kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
  kappa.chain[i,] <- kappa
  
  ### CONSTRAINTS
  level <- mean(kappa.chain[i,])
  incline <- 1/sum(beta.chain[i-1,])
  #incline <- (kappa.chain[i,1] - kappa.chain[i,t]) / (t - 1)
  
  #kappa.chain[i,] <- (kappa.chain[i, ] - level) / incline
  
  level2 <- mean(theta.chain[i-1,])
  incline2 <- 1/sum(gamma.chain[i-1,])
  
  theta.chain[i-1,] <- (theta.chain[i-1,] - level2) / incline2
  
  alpha.chain[i-1,] <- alpha.chain[i-1,] + (beta.chain[i-1,] * level) + (gamma.chain[i-1,] * level2)
  beta.chain[i-1,] <- beta.chain[i-1,] * incline
  gamma.chain[i-1,] <- gamma.chain[i-1,] * incline2
  
  
  
  ### SIGMA_E
  # A <- matrix(kronecker(rep(1, n), alpha.chain[i-1, ]))
  # B <- matrix(kronecker(rep(1, n), beta.chain[i-1, ]))
  # C <- matrix(kronecker(theta.chain[i-1, ], gamma.chain[i-1, ]))
  # 
  # aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i, ]) - C[, rep(1, t)])^2 )
  # sigma_e.chain[i] <- invgamma::rinvgamma(1, t*n*ages/2, rate = aux/2) 
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### ALPHA, BETA AND GAMMA
  for (j in 1:ages) {
    Y.aux <- c( y_obs[j, ] )
    for(k in 2:n){
      Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
    }
    
    A <- rep(alpha.chain[i-1, j], n*t)
    B <- rep(beta.chain[i-1, j] * kappa.chain[i,], n)
    C <- rep(gamma.chain[i-1, j] * theta.chain[i-1,], t)
    
    aux <- sum( (Y.aux - A - B - C)^2 )
    sigma_e.chain[i, j] <- invgamma::rinvgamma(1, t*n/2, rate = aux/2)
    
    
    X <- cbind(1, rep(kappa.chain[i, ], n), rep(theta.chain[i-1, ], each = t))
    aux.reg <- chol2inv(chol(t(X) %*% X))
    mean.reg <- aux.reg %*% t(X) %*% matrix(Y.aux)
    var.reg <- sigma_e.chain[i, j] * aux.reg
    
    tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
    
    alpha.chain[i, j] <- tmp[1]
    beta.chain[i, j] <- tmp[2]
    gamma.chain[i, j] <- tmp[3]
  }
  
  ### BETA & GAMMA CONSTRAINT
  # beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])
  # gamma.chain[i, ] <- gamma.chain[i, ]/sum(gamma.chain[i, ])
  
  ### SIGMA_W
  aux <- sum( (kappa.chain[i, 2:t] - nu.chain[i-1] - kappa.chain[i, 1:(t-1)])^2 ) ## years
  sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### NU
  aux <- sigma_w.chain[i]/(t-1)
  aux.2 <- (kappa.chain[i, t] - kappa.chain[i, 1])/(t-1)
  nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
  #nu.chain[i] <- -0.95
  
  ### SIGMA_T
  aux <- Mi - Wi
  aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
  sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
  
  ### THETA_q
  for(j in 1:n){
    mu0 = sum( (Wi[j, -j]/Mi[j,j])*theta.chain[i-1, -j] )
    sigma0 = sigma_t.chain[i]/Mi[j,j]
    
    
    y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta.chain[i,]) %*% kappa.chain[i, ]) - matrix(alpha.chain[i,],
                                                                                                    nrow = ages,
                                                                                                    ncol = t) )
    
    #aux = ( sigma0*(t*sum(gamma.chain[i,]^2)) + sigma_e.chain[i] )
    aux = ( sigma0*t*sum( (gamma.chain[i,]^2)/sigma_e.chain[i,] ) + 1 )
    #aux.2 = y.hat*matrix(gamma.chain[i,], nrow = ages, ncol = t)
    aux.2 = y.hat*matrix((gamma.chain[i,]/sigma_e.chain[i,]), nrow = ages, ncol = t)
    #mu1 = ( sigma0*sum(aux.2) + sigma_e.chain[i]*mu0 )/aux
    mu1 = ( sigma0*sum(aux.2) + mu0 )/aux
    #sigma1 = (sigma_e.chain[i]*sigma0)/aux
    sigma1 = sigma0/aux
    
    theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
  }
  
  
}
## 2m 10000 iterations

### chain treatment
bn = 10000
thin = 5
kappa.est <- apply(fit$kappa[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

alpha.est <- apply(fit$alpha[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
beta.est <- apply(fit$beta[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
gamma.est <- apply(fit$gamma[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
theta.est <- apply(fit$theta[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

sigma_w.est <- quantile(fit$sigma_w[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_t.est <- quantile(fit$sigma_t[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_e.est <- apply(fit$sigma_e[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975))
nu.est <- quantile(fit$nu[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
#lambda.est <- quantile(lambda.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))

#### age-time parameters plot
x11()
par(mfrow = c(2,2))
plot(1:t, kappa)
lines(1:t, kappa.est[2,])
lines(1:t, kappa.est[1,], lty = 2, col = "blue")
lines(1:t, kappa.est[3,], lty = 2, col = "blue")

plot(1:ages, alpha)
lines(1:ages, alpha.est[2,])
lines(1:ages, alpha.est[1,], lty = 2, col = "blue")
lines(1:ages, alpha.est[3,], lty = 2, col = "blue")

plot(1:ages, beta)
lines(1:ages, beta.est[2,])
lines(1:ages, beta.est[1,], lty = 2, col = "blue")
lines(1:ages, beta.est[3,], lty = 2, col = "blue")

plot(1:ages, gamma)
lines(1:ages, gamma.est[2,])
lines(1:ages, gamma.est[1,], lty = 2, col = "blue")
lines(1:ages, gamma.est[3,], lty = 2, col = "blue")

plot(1:n, theta)
lines(1:n, theta.est[2,])
lines(1:n, theta.est[1,], lty = 2, col = "blue")
lines(1:n, theta.est[3,], lty = 2, col = "blue")

# plot.ts(kappa.chain[bn:it ,1])
# plot.ts(kappa.chain[bn:it ,4])
# plot.ts(kappa.chain[bn:it ,7])
# plot.ts(kappa.chain[bn:it ,10])
# 
# acf(kappa.chain[bn:it ,1])
# acf(kappa.chain[bn:it ,4])
# acf(kappa.chain[bn:it ,7])
# acf(kappa.chain[bn:it ,10])

plot.ts(alpha.chain[seq(bn, it, by = thin),1])
plot.ts(fit$alpha[bn:it ,4])
plot.ts(fit$alpha[bn:it ,8])
plot.ts(fit$alpha[bn:it ,12])

acf(alpha.chain[bn:it ,1])
acf(alpha.chain[bn:it ,4])
acf(alpha.chain[bn:it ,8])
acf(alpha.chain[bn:it ,12])

plot.ts(fit$beta[bn:it, 1])
plot.ts(beta.chain[bn:it, 4])
plot.ts(beta.chain[bn:it, 8])
plot.ts(beta.chain[bn:it, 12])

acf(beta.chain[bn:it, 1])
acf(beta.chain[bn:it, 4])
acf(beta.chain[bn:it, 8])
acf(beta.chain[bn:it, 12])

plot.ts(fit$gamma[bn:it, 1])
plot.ts(gamma.chain[bn:it, 4])
plot.ts(gamma.chain[bn:it, 8])
plot.ts(gamma.chain[bn:it, 12])

acf(gamma.chain[bn:it, 1])
acf(gamma.chain[bn:it, 4])
acf(gamma.chain[bn:it, 8])
acf(gamma.chain[bn:it, 12])

plot.ts(fit$theta[bn:it, 1])
plot.ts(theta.chain[bn:it, 4])
plot.ts(theta.chain[bn:it, 8])
plot.ts(theta.chain[bn:it, 12])

acf(theta.chain[bn:it, 1])
acf(theta.chain[bn:it, 5])
acf(theta.chain[bn:it, 10])
acf(theta.chain[bn:it, 15])
graphics.off()


nu.est
sigma_w.est
sigma_e.est
sigma_t.est
theta.est
#lambda.est

#par(mfrow = c(1, 3))
plot.ts(nu.chain[seq(bn, it, by = thin)])
abline(h = nu, col = 2)
plot.ts(sigma_w.chain[seq(bn, it, by = thin)])
abline(h = sigma_w, col = 2)

plot(sigma_e, ylim = c(0, 1))
lines(sigma_e.est[2,])
lines(sigma_e.est[1,], col = "blue", lty = 2)
lines(sigma_e.est[3,], col = "blue", lty = 2)
### ok

plot.ts(sigma_t.chain[seq(bn, it, by = thin)])
abline(h = sigma_t, col = 2)

par(mfrow = c(1,3))
plot.ts(theta.chain[seq(bn, it, by = thin), 1])
abline(h = theta[1], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 5])
abline(h = theta[5], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 10])
abline(h = theta[10], col = 2)
# plot.ts(lambda.chain[seq(bn, it, by = thin)])
# abline(h = lambda, col = 2)

acf(nu.chain[bn:it])
acf(sigma_w.chain[bn:it])
acf(sigma_e.chain[bn:it])
acf(sigma_t.chain[bn:it])
graphics.off()

fit <- list(kappa = kappa.chain,
            alpha = alpha.chain,
            beta = beta.chain,
            gamma = gamma.chain,
            theta = theta.chain,
            sigma_w = sigma_w.chain,
            sigma_e = sigma_e.chain,
            sigma_t = sigma_t.chain,
            #lambda = lambda.chain,
            nu = nu.chain)

saveRDS(fit, "spatialLCfit_agesigmae.RDS")
fit <- readRDS("spatialLCfit_agesigmae.RDS")
### estimation (n-steps) #age-time-spatial parameters (gamma + theta constraint)----
it = 20000 #iterations
kappa.chain <- matrix(NA, nrow = it, ncol = t)

alpha.chain <- matrix(NA, nrow = it, ncol = ages)
beta.chain <- matrix(NA, nrow = it, ncol = ages)
gamma.chain <- matrix(NA, nrow = it, ncol = ages)

theta.chain <- matrix(NA, nrow = it, ncol = n)

sigma_w.chain <- rep(NA, it)
sigma_e.chain <- rep(NA, it)
sigma_t.chain <- rep(NA, it)
#lambda.chain <- rep(NA, it)
nu.chain <- rep(NA, it)

## starters
alpha.chain[1, ] <- runif(ages)
#alpha.chain[1, ] <- alpha.chain[1, ]/sum(alpha.chain[1, ])
beta.chain[1, ] <- runif(ages)
beta.chain[1, ] <- beta.chain[1, ]/sum(beta.chain[1, ])
gamma.chain[1, ] <- runif(ages)
gamma.chain[1, ] <- gamma.chain[1, ]/sum(gamma.chain[1, ])

theta.chain[1, ] <- theta    ### theta = 1 1 1 causou instabilidade

sigma_w.chain[1] <- 1
sigma_e.chain[1] <- 1
sigma_t.chain[1] <- 1

#lambda.min <- 1/min(labda.dom.i$values); lambda.max <- 1/max(labda.dom.i$values)
#lambda.chain[1] <- runif(1, lambda.min, lambda.max)
#lambda.chain[1] = lambda

nu.chain[1] <- runif(1, -1, 1)

## gibbs
pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = it, clear = FALSE, width = 60)
pb$tick()
prop = 0

set.seed(16)
for(i in 2:it){
  
  pb$tick()
  ### kappa
  fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                 sigma_e = sigma_e.chain[i-1],
                 sigma_w = sigma_w.chain[i-1],
                 nu = nu.chain[i-1],
                 alpha = alpha.chain[i-1, ],
                 beta = beta.chain[i-1, ],
                 gamma = gamma.chain[i-1, ],
                 theta = theta.chain[i-1, ])
  
  kappa.chain[i,] <- fit$kappa.bs
  ### kappa constraint
  kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
  #kappa.chain[i,] <- kappa
  
  ### sigma_e
  A <- matrix(kronecker(rep(1, n), alpha.chain[i-1, ]))
  B <- matrix(kronecker(rep(1, n), beta.chain[i-1, ]))
  C <- matrix(kronecker(theta.chain[i-1, ], gamma.chain[i-1, ]))
  #C <- matrix(kronecker(theta, gamma.chain[i-1, ]))
  
  aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i,]) - C[, rep(1, t)])^2 )
  sigma_e.chain[i] <- invgamma::rinvgamma(1, t*n*ages/2, rate = aux/2) ##gerar fixando idade?
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### alpha, beta and gamma 
  for (j in 1:ages) {
    Y.aux <- c( y_obs[j, ] )
    for(k in 2:n){
      Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
    }
    
    X <- cbind(1, rep(kappa.chain[i,], n), rep(theta.chain[i-1, ], each = t))
    #X <- cbind(1, rep(kappa.chain[i,], n), rep(theta, each = t))
    aux.reg <- chol2inv(chol(t(X) %*% X))
    mean.reg <- aux.reg %*% t(X) %*% matrix(Y.aux)
    var.reg <- sigma_e.chain[i] * aux.reg
    
    tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
    
    alpha.chain[i, j] <- tmp[1]
    beta.chain[i, j] <- tmp[2]
    gamma.chain[i, j] <- tmp[3]
  }
  
  ### alpha, beta & gamma constraint
  #alpha.chain[i, ] <- alpha.chain[i, ]/sum(alpha.chain[i, ])
  beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])
  gamma.chain[i, ] <- gamma.chain[i, ]/sum(gamma.chain[i, ])
  
  ### sigma_w
  aux <- sum( (kappa.chain[i, 2:t] - nu.chain[i-1] - kappa.chain[i, 1:(t-1)])^2 ) ## years
  sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### nu
  aux <- sigma_w.chain[i]/(t-1)
  aux.2 <- (kappa.chain[i, t] - kappa.chain[i, 1])/(t-1)
  nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
  #nu.chain[i] <- nu
  
  ### sigma_t
  #aux <- Mi - lambda.chain[i-1]*Wi
  aux <- Mi - Wi
  aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
  #aux.2 <- theta %*% aux %*% matrix(theta)
  sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
  ### um pouco esquisito valores mt altos surgirem da estimacao
  
  ### theta_q 
  for(j in 1:n){
    mu0 = sum( (Wi[j, -j]/Mi[j,j])*theta.chain[i-1, -j] )
    ### é a de cima menos lambda
    #mu0 = sum( Wi[j, -j] * (theta.chain[i-1, j] - theta.chain[i-1, -j])^2 ) 
    sigma0 = sigma_t.chain[i]/Mi[j,j]
    #sigma0 = sigma_t.chain[i]
    
    y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta.chain[i,]) %*% kappa.chain[i,]) - matrix(alpha.chain[i,],
                                                                                                   nrow = ages,
                                                                                                   ncol = t) )
    
    aux = ( sigma0*(t*sum(gamma.chain[i,]^2)) + sigma_e.chain[i] )
    aux.2 = y.hat*matrix(gamma.chain[i,], nrow = ages, ncol = t)
    mu1 = ( sigma0*sum(aux.2) + sigma_e.chain[i]*mu0 )/aux
    sigma1 = (sigma_e.chain[i]*sigma0)/aux
    
    theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
  }
  
  theta.chain[i,] <- theta.chain[i,] - mean(theta.chain[i,])
  
  
  
}
## 2m 10000 iterations

### chain treatment
bn = 10000
thin = 5
kappa.est <- apply(kappa.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

alpha.est <- apply(alpha.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
beta.est <- apply(beta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
gamma.est <- apply(gamma.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
theta.est <- apply(theta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

sigma_w.est <- quantile(sigma_w.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_t.est <- quantile(sigma_t.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_e.est <- quantile(sigma_e.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
nu.est <- quantile(nu.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
#lambda.est <- quantile(lambda.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))

#### age-time parameters plot
x11()
par(mfrow = c(2,2))
plot(1:t, kappa)
lines(1:t, kappa.est[2,])
lines(1:t, kappa.est[1,], lty = 2, col = "blue")
lines(1:t, kappa.est[3,], lty = 2, col = "blue")

plot(1:ages, alpha)
lines(1:ages, alpha.est[2,])
lines(1:ages, alpha.est[1,], lty = 2, col = "blue")
lines(1:ages, alpha.est[3,], lty = 2, col = "blue")

plot(1:ages, beta)
lines(1:ages, beta.est[2,])
lines(1:ages, beta.est[1,], lty = 2, col = "blue")
lines(1:ages, beta.est[3,], lty = 2, col = "blue")

plot(1:ages, gamma)
lines(1:ages, gamma.est[2,])
lines(1:ages, gamma.est[1,], lty = 2, col = "blue")
lines(1:ages, gamma.est[3,], lty = 2, col = "blue")

plot(1:n, theta)
lines(1:n, theta.est[2,])
lines(1:n, theta.est[1,], lty = 2, col = "blue")
lines(1:n, theta.est[3,], lty = 2, col = "blue")

# plot.ts(kappa.chain[bn:it ,1])
# plot.ts(kappa.chain[bn:it ,4])
# plot.ts(kappa.chain[bn:it ,7])
# plot.ts(kappa.chain[bn:it ,10])
# 
# acf(kappa.chain[bn:it ,1])
# acf(kappa.chain[bn:it ,4])
# acf(kappa.chain[bn:it ,7])
# acf(kappa.chain[bn:it ,10])

plot.ts(alpha.chain[seq(bn, it, by = thin),1])
plot.ts(alpha.chain[bn:it ,4])
plot.ts(alpha.chain[bn:it ,8])
plot.ts(alpha.chain[bn:it ,12])

acf(alpha.chain[bn:it ,1])
acf(alpha.chain[bn:it ,4])
acf(alpha.chain[bn:it ,8])
acf(alpha.chain[bn:it ,12])

plot.ts(beta.chain[bn:it, 1])
plot.ts(beta.chain[bn:it, 4])
plot.ts(beta.chain[bn:it, 8])
plot.ts(beta.chain[bn:it, 12])

acf(beta.chain[bn:it, 1])
acf(beta.chain[bn:it, 4])
acf(beta.chain[bn:it, 8])
acf(beta.chain[bn:it, 12])

plot.ts(gamma.chain[bn:it, 1])
plot.ts(gamma.chain[bn:it, 4])
plot.ts(gamma.chain[bn:it, 8])
plot.ts(gamma.chain[bn:it, 12])

acf(gamma.chain[bn:it, 1])
acf(gamma.chain[bn:it, 4])
acf(gamma.chain[bn:it, 8])
acf(gamma.chain[bn:it, 12])

plot.ts(theta.chain[bn:it, 1])
plot.ts(theta.chain[bn:it, 4])
plot.ts(theta.chain[bn:it, 8])
plot.ts(theta.chain[bn:it, 12])

acf(theta.chain[bn:it, 1])
acf(theta.chain[bn:it, 5])
acf(theta.chain[bn:it, 10])
acf(theta.chain[bn:it, 15])
graphics.off()


nu.est
sigma_w.est
sigma_e.est
sigma_t.est
theta.est
#lambda.est

#par(mfrow = c(1, 3))
plot.ts(nu.chain[seq(bn, it, by = thin)], ylim = c(-0.55, -0.7))
abline(h = nu, col = 2)
plot.ts(sigma_w.chain[seq(bn, it, by = thin)])
abline(h = sigma_w, col = 2)
plot.ts(sigma_e.chain[seq(bn, it, by = thin)])
abline(h = sigma_e, col = 2)
plot.ts(sigma_t.chain[seq(bn, it, by = thin)])
abline(h = sigma_t, col = 2)

par(mfrow = c(1,3))
plot.ts(theta.chain[seq(bn, it, by = thin), 1])
abline(h = theta[1], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 5])
abline(h = theta[5], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 10])
abline(h = theta[10], col = 2)
# plot.ts(lambda.chain[seq(bn, it, by = thin)])
# abline(h = lambda, col = 2)

acf(nu.chain[bn:it])
acf(sigma_w.chain[bn:it])
acf(sigma_e.chain[bn:it])
acf(sigma_t.chain[bn:it])
graphics.off()

fit <- list(kappa = kappa.chain,
            alpha = alpha.chain,
            beta = beta.chain,
            gamma = gamma.chain,
            theta = theta.chain,
            sigma_w = sigma_w.chain,
            sigma_e = sigma_e.chain,
            sigma_t = sigma_t.chain,
            #lambda = lambda.chain,
            nu = nu.chain)

saveRDS(fit, "spatialLCfit_nolambda_nok0.RDS")
### sigma_t estranho
### estimation (n-steps) #age-time-spatial parameters (k0, gamma + theta constraint) ----
it = 20000 #iterations
kappa.chain <- matrix(NA, nrow = it, ncol = t+1)

alpha.chain <- matrix(NA, nrow = it, ncol = ages)
beta.chain <- matrix(NA, nrow = it, ncol = ages)
gamma.chain <- matrix(NA, nrow = it, ncol = ages)

theta.chain <- matrix(NA, nrow = it, ncol = n)

sigma_w.chain <- rep(NA, it)
sigma_e.chain <- rep(NA, it)
sigma_t.chain <- rep(NA, it)
#lambda.chain <- rep(NA, it)
nu.chain <- rep(NA, it)

## starters
alpha.chain[1, ] <- runif(ages)
#alpha.chain[1, ] <- alpha.chain[1, ]/sum(alpha.chain[1, ])
beta.chain[1, ] <- runif(ages)
beta.chain[1, ] <- beta.chain[1, ]/sum(beta.chain[1, ])
gamma.chain[1, ] <- runif(ages)
gamma.chain[1, ] <- gamma.chain[1, ]/sum(gamma.chain[1, ])

theta.chain[1, ] <- theta    ### theta = 1 1 1 causou instabilidade

sigma_w.chain[1] <- 1
sigma_e.chain[1] <- 1
sigma_t.chain[1] <- 1

#lambda.min <- 1/min(labda.dom.i$values); lambda.max <- 1/max(labda.dom.i$values)
#lambda.chain[1] <- runif(1, lambda.min, lambda.max)
#lambda.chain[1] = lambda

nu.chain[1] <- runif(1, -1, 1)

## gibbs
pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = it, clear = FALSE, width = 60)
pb$tick()
prop = 0

set.seed(16)
for(i in 2:it){
  
  pb$tick()
  ### kappa
  fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                 sigma_e = sigma_e.chain[i-1],
                 sigma_w = sigma_w.chain[i-1],
                 nu = nu.chain[i-1],
                 alpha = alpha.chain[i-1, ],
                 beta = beta.chain[i-1, ],
                 gamma = gamma.chain[i-1, ],
                 theta = theta.chain[i-1, ])

  kappa.chain[i, 2:(t+1)] <- fit$kappa.bs
  kappa.chain[i, 1] <- fit$kappa0.bs
  ### kappa constraint
  kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
  
  ### sigma_e
  A <- matrix(kronecker(rep(1, n), alpha.chain[i-1, ]))
  B <- matrix(kronecker(rep(1, n), beta.chain[i-1, ]))
  C <- matrix(kronecker(theta.chain[i-1, ], gamma.chain[i-1, ]))
  #C <- matrix(kronecker(theta, gamma.chain[i-1, ]))
  
  aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i, 2:(t+1)]) - C[, rep(1, t)])^2 )
  sigma_e.chain[i] <- invgamma::rinvgamma(1, t*n*ages/2, rate = aux/2) ##gerar fixando idade?
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### alpha, beta and gamma 
  for (j in 1:ages) {
    Y.aux <- c( y_obs[j, ] )
    for(k in 2:n){
      Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
    }
    
    X <- cbind(1, rep(kappa.chain[i, 2:(t+1)], n), rep(theta.chain[i-1, ], each = t))
    #X <- cbind(1, rep(kappa.chain[i,], n), rep(theta, each = t))
    aux.reg <- chol2inv(chol(t(X) %*% X))
    mean.reg <- aux.reg %*% t(X) %*% matrix(Y.aux)
    var.reg <- sigma_e.chain[i] * aux.reg
    
    tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
    
    alpha.chain[i, j] <- tmp[1]
    beta.chain[i, j] <- tmp[2]
    gamma.chain[i, j] <- tmp[3]
  }
  
  ### alpha, beta & gamma constraint
  #alpha.chain[i, ] <- alpha.chain[i, ]/sum(alpha.chain[i, ])
  beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])
  gamma.chain[i, ] <- gamma.chain[i, ]/sum(gamma.chain[i, ])
  
  ### sigma_w
  aux <- sum( (kappa.chain[i, 2:(t+1)] - nu.chain[i-1] - kappa.chain[i, 1:t])^2 ) ## years
  sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### nu
  aux <- sigma_w.chain[i]/(t-1)
  aux.2 <- (kappa.chain[i, t+1] - kappa.chain[i, 1])/t
  nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
  #nu.chain[i] <- nu
  
  ### sigma_t
  #aux <- Mi - lambda.chain[i-1]*Wi
  aux <- Mi - Wi
  aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
  #aux.2 <- theta %*% aux %*% matrix(theta)
  sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
  ### um pouco esquisito valores mt altos surgirem da estimacao
  
  ### theta_q 
  for(j in 1:n){
    mu0 = sum( (Wi[j, -j]/Mi[j,j])*theta.chain[i-1, -j] )
    ### é a de cima menos lambda
    #mu0 = sum( Wi[j, -j] * (theta.chain[i-1, j] - theta.chain[i-1, -j])^2 ) 
    sigma0 = sigma_t.chain[i]/Mi[j,j]
    #sigma0 = sigma_t.chain[i]
    
    y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta.chain[i,]) %*% kappa.chain[i, 2:(t+1)]) - matrix(alpha.chain[i,],
                                                                                                   nrow = ages,
                                                                                                   ncol = t) )
    
    aux = ( sigma0*(t*sum(gamma.chain[i,]^2)) + sigma_e.chain[i] )
    aux.2 = y.hat*matrix(gamma.chain[i,], nrow = ages, ncol = t)
    mu1 = ( sigma0*sum(aux.2) + sigma_e.chain[i]*mu0 )/aux
    sigma1 = (sigma_e.chain[i]*sigma0)/aux
    
    theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
  }
  
  theta.chain[i,] <- theta.chain[i,] - mean(theta.chain[i,])
  
  
  
}
## 2m 10000 iterations

### chain treatment
bn = 10000
thin = 5
kappa.est <- apply(kappa.chain[seq(bn, it, by = thin), 2:41], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

alpha.est <- apply(alpha.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
beta.est <- apply(beta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
gamma.est <- apply(gamma.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
theta.est <- apply(theta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

sigma_w.est <- quantile(sigma_w.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_t.est <- quantile(sigma_t.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_e.est <- quantile(sigma_e.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
nu.est <- quantile(nu.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
#lambda.est <- quantile(lambda.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))

#### age-time parameters plot
x11()
par(mfrow = c(2,2))
plot(1:t, kappa)
lines(1:t, kappa.est[2, ])
lines(1:t, kappa.est[1, ], lty = 2, col = "blue")
lines(1:t, kappa.est[3, ], lty = 2, col = "blue")

plot(1:ages, alpha)
lines(1:ages, alpha.est[2,])
lines(1:ages, alpha.est[1,], lty = 2, col = "blue")
lines(1:ages, alpha.est[3,], lty = 2, col = "blue")

plot(1:ages, beta)
lines(1:ages, beta.est[2,])
lines(1:ages, beta.est[1,], lty = 2, col = "blue")
lines(1:ages, beta.est[3,], lty = 2, col = "blue")

plot(1:ages, gamma)
lines(1:ages, gamma.est[2,])
lines(1:ages, gamma.est[1,], lty = 2, col = "blue")
lines(1:ages, gamma.est[3,], lty = 2, col = "blue")

plot(1:n, theta)
lines(1:n, theta.est[2,])
lines(1:n, theta.est[1,], lty = 2, col = "blue")
lines(1:n, theta.est[3,], lty = 2, col = "blue")

plot.ts(kappa.chain[bn:it ,1])
plot.ts(kappa.chain[bn:it ,4])
plot.ts(kappa.chain[bn:it ,7])
plot.ts(kappa.chain[bn:it ,10])
# 
# acf(kappa.chain[bn:it ,1])
# acf(kappa.chain[bn:it ,4])
# acf(kappa.chain[bn:it ,7])
# acf(kappa.chain[bn:it ,10])

plot.ts(alpha.chain[seq(bn, it, by = thin),1])
plot.ts(alpha.chain[bn:it ,4])
plot.ts(alpha.chain[bn:it ,8])
plot.ts(alpha.chain[bn:it ,12])

acf(alpha.chain[bn:it ,1])
acf(alpha.chain[bn:it ,4])
acf(alpha.chain[bn:it ,8])
acf(alpha.chain[bn:it ,12])

plot.ts(beta.chain[bn:it, 1])
plot.ts(beta.chain[bn:it, 4])
plot.ts(beta.chain[bn:it, 8])
plot.ts(beta.chain[bn:it, 12])

acf(beta.chain[bn:it, 1])
acf(beta.chain[bn:it, 4])
acf(beta.chain[bn:it, 8])
acf(beta.chain[bn:it, 12])

plot.ts(gamma.chain[bn:it, 1])
plot.ts(gamma.chain[bn:it, 4])
plot.ts(gamma.chain[bn:it, 8])
plot.ts(gamma.chain[bn:it, 12])

acf(gamma.chain[bn:it, 1])
acf(gamma.chain[bn:it, 4])
acf(gamma.chain[bn:it, 8])
acf(gamma.chain[bn:it, 12])

plot.ts(theta.chain[bn:it, 1])
plot.ts(theta.chain[bn:it, 4])
plot.ts(theta.chain[bn:it, 8])
plot.ts(theta.chain[bn:it, 12])

acf(theta.chain[bn:it, 1])
acf(theta.chain[bn:it, 5])
acf(theta.chain[bn:it, 10])
acf(theta.chain[bn:it, 15])
graphics.off()


nu.est
sigma_w.est
sigma_e.est
sigma_t.est
theta.est
#lambda.est

#par(mfrow = c(1, 3))
plot.ts(nu.chain[seq(bn, it, by = thin)], ylim = c(-0.55, -0.7))
abline(h = nu, col = 2)
plot.ts(sigma_w.chain[seq(bn, it, by = thin)])
abline(h = sigma_w, col = 2)
plot.ts(sigma_e.chain[seq(bn, it, by = thin)])
abline(h = sigma_e, col = 2)
plot.ts(sigma_t.chain[seq(bn, it, by = thin)])
abline(h = sigma_t, col = 2)

par(mfrow = c(1,3))
plot.ts(theta.chain[seq(bn, it, by = thin), 1])
abline(h = theta[1], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 5])
abline(h = theta[5], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 10])
abline(h = theta[10], col = 2)
# plot.ts(lambda.chain[seq(bn, it, by = thin)])
# abline(h = lambda, col = 2)

acf(nu.chain[bn:it])
acf(sigma_w.chain[bn:it])
acf(sigma_e.chain[bn:it])
acf(sigma_t.chain[bn:it])
graphics.off()

fit <- list(kappa = kappa.chain,
            alpha = alpha.chain,
            beta = beta.chain,
            gamma = gamma.chain,
            theta = theta.chain,
            sigma_w = sigma_w.chain,
            sigma_e = sigma_e.chain,
            sigma_t = sigma_t.chain,
            #lambda = lambda.chain,
            nu = nu.chain)

saveRDS(fit, "spatialLCfit_nolambda.RDS")


##### MISSING COM DADO SIMULADO E FUNCAO SBLC_MISSING ----
# source("sblc_fun.R")
# source("sffbs_fun_k0.R")
# fit1 <- sblc(y_obs, ages, t, n, m0, C0, Mi, Wi, it = 2000, bn = 1000, thin = 1)
# median(fit1$nu.chain)
# nu

source("sblc_fun.R")
source("sffbs_fun.R")
fit2 <- sblc_2(y_obs, ages, t, n, m0, C0, Mi, Wi, it = 2000, bn = 1000, thin = 1)
median(fit2$nu.chain)
plot.ts(fit2$nu.chain)
plot.ts(fit2$kappa.chain[,1])
plot.ts(fit2$sigma_w.chain)
median(fit2$sigma_w.chain) ### alto
median(fit2$sigma_t.chain)  ###baixo
median(fit2$sigma_e.chain) ##bem estimado


# fit3 <- sblc_3(y_obs, ages, t, n, m0, C0, Mi, Wi, it = 2000, bn = 1000, thin = 1)
# median(fit3$nu.chain)
# plot.ts(fit3$nu.chain)
# plot.ts(fit3$kappa.chain[,1])
# plot.ts(fit3$sigma_w.chain)
# nu    ###ruim mudar de lugar

source("fitted_sblc.R")
fitted2 <- fitted_sblc(fit2)
aux.qx2 <- 1 - exp(-exp(fitted2$mean))
y_qx <- 1 - exp(-exp(y_obs))
res <- ( (y_qx - aux.qx2)^2)/aux.qx2

fitted3 <- fitted_sblc(fit3)
aux.qx3 <- 1 - exp(-exp(fitted3$mean))
res2 <- ( (y_qx - aux.qx3)^2)/aux.qx3
### fit 3 com maior residuo 

###bons valores, sigma_w e kappa nao convergiram mt bem, k0 teve desempenho pior

### MISSING VS INPUT MISSING
aux <- cbind(sample(1:195, 100, replace=T),
             sample(1:40, 100, replace = T))
y_obs[aux] <- NA

fit3 <- sblc_missing(y_obs, ages, t, n, m0, C0, Mi, Wi, it = 3000, bn = 1500, thin = 1)
median(fit3$nu.chain)
nu
