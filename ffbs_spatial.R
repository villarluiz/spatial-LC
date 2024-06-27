

ff_sp <- function(y_obs, ages, t, n, m0, C0, sigma_e, sigma_w, nu, alpha, beta, gamma, theta){
  ### aux ----
  A <- kronecker(rep(1, n), alpha)
  B <- kronecker(rep(1, n), beta)
  C <- kronecker(theta, gamma)
  
  ### estimation ----
  #### FF
  at <- Rt <- mt <- Ct <- rep(NA, t)
  ft <- matrix(NA, ncol = t, nrow = ages*n)
  Qt <- array(NA, dim=c(n*ages, n*ages, t))
  
  at[1] <- nu + m0
  Rt[1] <- C0 + sigma_w ##phiw
  
  ft[,1] <- A + (B * at[1]) + C
  Qt[,,1] <- Rt[1] * B %*% t(B) + diag(n*ages)*sigma_e #phie
  
  mt[1] <- at[1] + Rt[1] %*% t(B) %*% solve(Qt[,,1]) %*% (y_obs[,1] - ft[,1])
  Ct[1] <- Rt[1] - Rt[1] %*% t(B) %*% solve(Qt[,,1]) %*% B %*% Rt[1]
  for(k in 2:t){
    at[k] <- nu + m0
    Rt[k] <- C0 + sigma_w ##phiw
    
    ft[,k] <- A + (B * at[k]) + C
    Qt[,,k] <- Rt[k] * B %*% t(B) + diag(n*ages)*sigma_e #phie
    
    mt[k] <- at[k] + Rt[k] %*% t(B) %*% solve(Qt[,,k]) %*% (y_obs[,k] - ft[,k])
    Ct[k] <- Rt[k] - Rt[k] %*% t(B) %*% solve(Qt[,,k]) %*% B %*% Rt[k]
  }
  return(list(mt = mt,
              Ct = Ct,
              Rt = Rt,
              at = at,
              ft = ft,
              Qt = Qt))
}

bs_sp <- function(mt, Ct, at, Rt){
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
  return(list(kappa.bs = kappa.bs,
              st = st,
              St = St))
}

ffbs_sp <- function(y_obs, ages, t, n, m0, C0, sigma_e, sigma_w, nu, alpha, beta, gamma, theta){
  aux.ff <- ff_sp(y_obs, ages, t, n, m0, C0, sigma_e, sigma_w, nu, alpha, beta, gamma, theta)
  aux.bs <- bs_sp(aux.ff$mt, aux.ff$Ct, aux.ff$at, aux.ff$Rt)
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
labda.dom <- eigen(W) #### contem 1, por ex.
lambda = 0.9
Q = M - lambda*W

### assumindo apenas 3 microreg pra simplificar:
ages <- 40  ## idades
t <- 20  ## tempo
n <- 15   ## regioes

Qi <- Q[1:n, 1:n]
Wi <- W[1:n, 1:n]
Mi <- diag(rowSums(Wi))
labda.dom.i <- eigen(Wi) 

sigma_t = 0.5
set.seed(14)
theta <- MASS::mvrnorm(n = 1, mu = rep(0, n), Sigma = MASS::ginv(Qi)*sigma_t) ### sigma_theta = 1
#theta <- theta - mean(theta)

### alpha, beta and kappa ----
## acho que tem constraint demais...
set.seed(14)
alpha <- runif(ages, -1, 1)
#alpha <- alpha/sum(alpha)
beta <- runif(ages, -1, 1)
beta <- beta/sum(beta)
gamma <- runif(ages, -1, 1)
#gamma <- gamma/sum(gamma)
kappa0 <- 1
kappa <- rep(NA, t)
nu <- -0.5

set.seed(14)
kappa[1] <- nu + kappa0 + rnorm(1, 0, sqrt(0.025)) 
for(k in 2:t){
  kappa[k] <- nu + kappa[k-1] + rnorm(1, 0, sqrt(0.025))
}
kappa <- kappa - mean(kappa)

m0 <- 0; C0 <- 100
sigma_e = 0.05
sigma_w = 0.025

### aux ----
A <- kronecker(rep(1, n), alpha)
B <- kronecker(rep(1, n), beta)
C <- kronecker(theta, gamma)
set.seed(14)
y_obs <- matrix(NA, nrow = ages*n, ncol = t)
for(k in 1:t){
  y_obs[,k] <- A + (B * kappa[k]) + C + rnorm(n*ages, 0, sd = sqrt(0.05))
}


### estimation (1-step) #kappa ----
fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
               sigma_e = sigma_e, sigma_w = sigma_w, nu = nu, alpha = alpha,
               beta = beta, gamma = gamma, theta = theta)

plot(1:t, kappa)
lines(1:t, fit$kappa.bs)


### estimation (n-steps) #age-time parameters (after filtering) ----
M = 10000 #iterations
kappa.chain <- matrix(NA, nrow = M, ncol = t)

alpha.chain <- matrix(NA, nrow = M, ncol = ages)
beta.chain <- matrix(NA, nrow = M, ncol = ages)
gamma.chain <- matrix(NA, nrow = M, ncol = ages)

sigma_w.chain <- rep(NA, M)
sigma_e.chain <- rep(NA, M)
nu.chain <- rep(NA, M)

## starters
alpha.chain[1, ] <- runif(ages)
#alpha.chain[1, ] <- alpha.chain[1, ]/sum(alpha.chain[1, ])
beta.chain[1, ] <- runif(ages)
beta.chain[1, ] <- beta.chain[1, ]/sum(beta.chain[1, ])
gamma.chain[1, ] <- runif(ages)
#gamma.chain[1, ] <- gamma.chain[1, ]/sum(gamma.chain[1, ])

sigma_w.chain[1] <- 1
sigma_e.chain[1] <- 1
nu.chain[1] <- runif(1, -1, 1)

## gibbs
start <- proc.time()[[3]]
set.seed(16)
for(i in 2:M){
  
  ### kappa
  fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                 sigma_e = sigma_e.chain[i-1],
                 sigma_w = sigma_w.chain[i-1],
                 nu = nu.chain[i-1],
                 alpha = alpha.chain[i-1, ],
                 beta = beta.chain[i-1, ],
                 gamma = gamma.chain[i-1, ],
                 theta = theta)
  
  kappa.chain[i,] <- fit$kappa.bs
  ### kappa constraint
  kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
  
  ### sigma_e
  A <- matrix(kronecker(rep(1, n), alpha.chain[i-1, ]))
  B <- matrix(kronecker(rep(1, n), beta.chain[i-1, ]))
  C <- matrix(kronecker(theta, gamma.chain[i-1, ]))
  
  aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i,]) - C[, rep(1, t)])^2 )
  sigma_e.chain[i] <- invgamma::rinvgamma(1, 0.01 + t*n*ages/2, rate = 0.01 + aux/2) ##gerar fixando idade?
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### alpha, beta and gamma  (alpha e beta estao se perdendo!)
  for (j in 1:ages) {
    Y.aux <- c( y_obs[j, ] )
    for(k in 2:n){
      Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
    }
    
    X <- cbind(1, rep(kappa.chain[i,], n), rep(theta, each = t))
    aux.reg <- chol2inv(chol(t(X) %*% X))
    mean.reg <- aux.reg %*% t(X) %*% Y.aux
    var.reg <- sigma_e.chain[i] * aux.reg
    
    tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
    
    alpha.chain[i, j] <- tmp[1]
    beta.chain[i, j] <- tmp[2]
    gamma.chain[i, j] <- tmp[3]
  }
  
  ### beta constraint
  beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])   ## sum(beta) = 1
  #alpha.chain[i, ] <- alpha.chain[i, ]/sum(alpha.chain[i, ]) ## sum(alpha) = 1
  #gamma.chain[i, ] <- gamma.chain[i, ]/sum(gamma.chain[i, ]) ## sum(gamma) = 1
  
  ### sigma_w
  aux <- sum( (kappa.chain[i, 2:t] - nu.chain[i-1] - kappa.chain[i, 1:(t-1)])^2 ) ## years
  sigma_w.chain[i] <- invgamma::rinvgamma(1, 0.01 + (t-1)/2, rate = 0.01 + aux/2)
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### nu
  aux <- sigma_w.chain[i]/t            #(t-1)
  aux.2 <- (kappa.chain[i, t] - kappa.chain[i, 1])/t  #(t-1)
  nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
  
}
end <- proc.time()[[3]]
end - start ## seconds

### chain treatment
bn = 5000
thin = 5
it = M
kappa.est <- apply(kappa.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

alpha.est <- apply(alpha.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
beta.est <- apply(beta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
gamma.est <- apply(gamma.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

sigma_w.est <- quantile(sigma_w.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_e.est <- quantile(sigma_e.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
nu.est <- quantile(nu.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))

#### age-time parameters plot
x11()
par(mfrow = c(2,2))
plot(1:10, kappa)
lines(1:10, kappa.est[2,])
lines(1:10, kappa.est[1,], lty = 2, col = "blue")
lines(1:10, kappa.est[3,], lty = 2, col = "blue")

plot(1:20, alpha)
lines(1:20, alpha.est[2,])
lines(1:20, alpha.est[1,], lty = 2, col = "blue")
lines(1:20, alpha.est[3,], lty = 2, col = "blue")

plot(1:20, beta)
lines(1:20, beta.est[2,])
lines(1:20, beta.est[1,], lty = 2, col = "blue")
lines(1:20, beta.est[3,], lty = 2, col = "blue")

plot(1:20, gamma)
lines(1:20, gamma.est[2,])
lines(1:20, gamma.est[1,], lty = 2, col = "blue")
lines(1:20, gamma.est[3,], lty = 2, col = "blue")

plot.ts(kappa.chain[seq(bn, it, by = thin), 1])
plot.ts(kappa.chain[seq(bn, it, by = thin) ,4])
plot.ts(kappa.chain[seq(bn, it, by = thin) ,7])
plot.ts(kappa.chain[seq(bn, it, by = thin) ,9])

acf(kappa.chain[bn:M ,1])
acf(kappa.chain[bn:M ,4])
acf(kappa.chain[bn:M ,7])
acf(kappa.chain[bn:M ,10])

plot.ts(alpha.chain[bn:M ,1])
plot.ts(alpha.chain[bn:M ,8])
plot.ts(alpha.chain[bn:M ,15])
plot.ts(alpha.chain[bn:M ,20])

acf(alpha.chain[bn:M ,1])
acf(alpha.chain[bn:M ,8])
acf(alpha.chain[bn:M ,15])
acf(alpha.chain[bn:M ,20])

plot.ts(beta.chain[bn:M, 1])
plot.ts(beta.chain[bn:M, 7])
plot.ts(beta.chain[bn:M, 16])
plot.ts(beta.chain[bn:M, 20])

acf(beta.chain[bn:M, 1])
acf(beta.chain[bn:M, 7])
acf(beta.chain[bn:M, 16])
acf(beta.chain[bn:M, 20])

plot.ts(gamma.chain[bn:M, 1])
plot.ts(gamma.chain[bn:M, 4])
plot.ts(gamma.chain[bn:M, 12])
plot.ts(gamma.chain[bn:M, 20])

acf(gamma.chain[bn:M, 1])
acf(gamma.chain[bn:M, 4])
acf(gamma.chain[bn:M, 12])
acf(gamma.chain[bn:M, 20])
graphics.off()


nu.est
sigma_w.est
sigma_e.est

par(mfrow = c(1, 3))
plot.ts(nu.chain[seq(bn, it, by = thin)])
abline(h = nu, col = 2)
plot.ts(sigma_w.chain[seq(bn, it, by = thin)])
abline(h = sigma_w, col = 2)
plot.ts(sigma_e.chain[seq(bn, it, by = thin)])
abline(h = sigma_e, col = 2)

acf(nu.chain[bn:M])
acf(sigma_w.chain[bn:M])
acf(sigma_e.chain[bn:M])

### estimation (n-steps) #error parameters (after filtering) ----
it = 10000 #iterations
kappa.chain <- matrix(NA, nrow = it, ncol = t)

sigma_t.chain <- rep(NA, it)
sigma_e.chain <- rep(NA, it)
sigma_w.chain <- rep(NA, it)

## starters
sigma_t.chain[1] <- 1
sigma_e.chain[1] <- 1
sigma_w.chain[1] <- 1

## gibbs
start <- proc.time()[[3]]
set.seed(16)
for(i in 2:it){
  
  ### kappa
  fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                 sigma_e = sigma_e.chain[i-1],
                 sigma_w = sigma_w.chain[i-1],
                 nu = nu,
                 alpha = alpha,
                 beta = beta,
                 gamma = gamma,
                 theta = theta)
  #theta = theta)
  
  kappa.chain[i,] <- fit$kappa.bs
  ### kappa constraint
  kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
  
  ### standardize?
  
  ### sigma_e
  A <- matrix(kronecker(rep(1, n), alpha))
  B <- matrix(kronecker(rep(1, n), beta))
  C <- matrix(kronecker(theta, gamma))
  
  aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i,]) - C[, rep(1, t)])^2 )
  sigma_e.chain[i] <- invgamma::rinvgamma(1, t*n*ages/2, rate = aux/2) ##gerar fixando idade?
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### alpha, beta and gamma 
  # for (j in 1:ages) {
  #   Y.aux <- c( y_obs[j, ] )
  #   for(k in 2:n){
  #     Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
  #   }
  #   
  #   X <- cbind(1, rep(kappa.chain[i,], n), rep(theta.chain[i-1, ], each = t))
  #   #X <- cbind(1, rep(kappa.chain[i,], n), rep(theta, each = t))
  #   aux.reg <- chol2inv(chol(t(X) %*% X))
  #   mean.reg <- aux.reg %*% t(X) %*% Y.aux
  #   var.reg <- sigma_e.chain[i] * aux.reg
  #   
  #   tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
  #   
  #   alpha.chain[i, j] <- tmp[1]
  #   beta.chain[i, j] <- tmp[2]
  #   gamma.chain[i, j] <- tmp[3]
  # }
  
  # ### beta constraint
  # beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])
  
  ### sigma_w
  aux <- sum( (kappa.chain[i, 2:t] - nu - kappa.chain[i, 1:(t-1)])^2 ) ## years
  sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  # ### nu
  # aux <- sigma_w.chain[i]/t            #(t-1)
  # aux.2 <- (kappa.chain[i, t] - kappa.chain[i, 1])/t  #(t-1)
  # nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
  
  ### sigma_t
  aux <- Mi - lambda*Wi
  aux.2 <- theta %*% aux %*% matrix(theta)
  sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
  ### um pouco esquisito valores mt altos surgirem da estimacao
  
  ### theta_q     ## rever contas, constraint nao deu jeito
  # for(j in 1:n){
  #   mu0 = sum((Wi[j, ]/diag(Mi)[j])*theta.chain[i-1, ])
  #   sigma0 = sigma_t.chain[i]/diag(Mi)[j]
  #   y.hat = (y_obs[1:ages + ages*(j-1), ] - (beta %*% t(kappa.chain[i,]) + alpha))
  #   
  #   aux = (sigma0*t*sum(gamma^2) + sigma_e)
  #   mu1 = (sigma0 * sum(y.hat*gamma) + sigma_e*mu0)/aux
  #   sigma1 = (sigma_e*sigma0)/aux
  #   
  #   theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
  # }
  # 
  # theta.chain[i,] <- theta.chain[i,] - mean(theta.chain[i,])
  
  # mu0 = lambda * rowSums(matrix(Wi/diag(Mi), nrow = n) %*% diag(theta.chain[i-1, ]))
  # sigma0 = sigma_t/diag(Mi)
  # y.hat = (y_obs - (beta.chain[i,] %*% t(kappa.chain[i,]) + alpha.chain[i,]))
  ### lambda
  
  
}
end <- proc.time()[[3]]
end - start ## seconds

### chain treatment
bn = 5000
thin = 5
kappa.est <- apply(kappa.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

sigma_t.est <- quantile(sigma_t.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975), na.rm = T)
sigma_e.est <- quantile(sigma_e.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975), na.rm = T)
sigma_w.est <- quantile(sigma_w.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975), na.rm = T)

#### age-time parameters plot
plot(1:10, kappa)
lines(1:10, kappa.est[2,])
lines(1:10, kappa.est[1,], lty = 2, col = "blue")
lines(1:10, kappa.est[3,], lty = 2, col = "blue")

par(mfrow = c(2,2))
plot.ts(kappa.chain[bn:it ,1])
plot.ts(kappa.chain[bn:it ,4])
plot.ts(kappa.chain[bn:it ,7])
plot.ts(kappa.chain[bn:it ,10])
graphics.off()


plot.ts(sigma_t.chain[seq(bn, it, by = thin)])
abline(h = sigma_t, col = 2)

plot.ts(sigma_w.chain[seq(bn, it, by = thin)])
abline(h = sigma_w, col = 2)

plot.ts(sigma_e.chain[seq(bn, it, by = thin)])
abline(h = sigma_e, col = 2)

### todos continuam estimando mal

### estimation (n-steps) #spatial parameters (no lambda, no sp error) (after filtering) ----
it = 10000 #iterations
kappa.chain <- matrix(NA, nrow = it, ncol = t)

theta.chain <- matrix(NA, nrow = it, ncol = n)

sigma_t.chain <- rep(NA, it)
#lambda.chain <- rep(NA, it)

## starters
theta.chain[1, ] <- theta    ### theta = 1 1 1 causou instabilidade

sigma_t.chain[1] <- 1

## gibbs
start <- proc.time()[[3]]
set.seed(16)
for(i in 2:it){
  
  ### kappa
  fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                 sigma_e = sigma_e,
                 sigma_w = sigma_w,
                 nu = nu,
                 alpha = alpha,
                 beta = beta,
                 gamma = gamma,
                 theta = theta.chain[i-1, ])
  #theta = theta)
  
  kappa.chain[i,] <- fit$kappa.bs
  ### kappa constraint
  kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
  
  ### standardize?
  
  # ### sigma_e
  # A <- matrix(kronecker(rep(1, n), alpha.chain[i-1, ]))
  # B <- matrix(kronecker(rep(1, n), beta.chain[i-1, ]))
  # C <- matrix(kronecker(theta.chain[i-1, ], gamma.chain[i-1, ]))
  # #C <- matrix(kronecker(theta, gamma.chain[i-1, ]))
  
  # aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i,]) - C[, rep(1, t)])^2 )
  # sigma_e.chain[i] <- invgamma::rinvgamma(1, t*n*ages/2, rate = aux/2) ##gerar fixando idade?
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### alpha, beta and gamma 
  # for (j in 1:ages) {
  #   Y.aux <- c( y_obs[j, ] )
  #   for(k in 2:n){
  #     Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
  #   }
  #   
  #   X <- cbind(1, rep(kappa.chain[i,], n), rep(theta.chain[i-1, ], each = t))
  #   #X <- cbind(1, rep(kappa.chain[i,], n), rep(theta, each = t))
  #   aux.reg <- chol2inv(chol(t(X) %*% X))
  #   mean.reg <- aux.reg %*% t(X) %*% Y.aux
  #   var.reg <- sigma_e.chain[i] * aux.reg
  #   
  #   tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
  #   
  #   alpha.chain[i, j] <- tmp[1]
  #   beta.chain[i, j] <- tmp[2]
  #   gamma.chain[i, j] <- tmp[3]
  # }
  
  # ### beta constraint
  # beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])
  
  # ### sigma_w
  # aux <- sum( (kappa.chain[i, 2:t] - nu.chain[i-1] - kappa.chain[i, 1:(t-1)])^2 ) ## years
  # sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
  # ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  # ### nu
  # aux <- sigma_w.chain[i]/t            #(t-1)
  # aux.2 <- (kappa.chain[i, t] - kappa.chain[i, 1])/t  #(t-1)
  # nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
  
  ### sigma_t
  aux <- Mi - lambda*Wi
  aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
  # #aux.2 <- theta %*% aux %*% matrix(theta)
  sigma_t.chain[i] <- invgamma::rinvgamma(1, 0.01 + (n-1)/2, rate = 0.01 + aux.2/2)
  ### um pouco esquisito valores mt altos surgirem da estimacao
  
  ### theta_q     ## rever contas, constraint nao deu jeito
  for(j in 1:n){
    mu0 = sum( (Wi[j, ]/Mi[j,j])*theta.chain[i-1, ] )
    sigma0 = sigma_t/Mi[j,j]
    
    y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta) %*% kappa.chain[i,]) - matrix(alpha,
                                                                                        nrow = ages,
                                                                                        ncol = t) )
    
    aux = ( sigma0*(t*sum(gamma^2)) + sigma_e )
    aux.2 = y.hat*matrix(gamma, nrow = ages, ncol = t)
    mu1 = ( sigma0*sum(aux.2) + sigma_e*mu0 )/aux
    sigma1 = (sigma_e*sigma0)/aux
    
    theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
  }
  
  theta.chain[i,] <- theta.chain[i,] - mean(theta.chain[i,])
  
  # mu0 = lambda * rowSums(matrix(Wi/diag(Mi), nrow = n) %*% diag(theta.chain[i-1, ]))
  # sigma0 = sigma_t/diag(Mi)
  # y.hat = (y_obs - (beta.chain[i,] %*% t(kappa.chain[i,]) + alpha.chain[i,]))
  ### lambda
  
  
}
end <- proc.time()[[3]]
end - start ## seconds

### chain treatment
bn = 5000
kappa.est <- apply(kappa.chain[bn:it, ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

sigma_t.est <- quantile(sigma_t.chain[bn:it], c(0.025, 0.5, 0.975), na.rm = T)
theta.est <- apply(theta.chain[bn:it,], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

#### age-time parameters plot
x11()
par(mfrow = c(2,2))
plot(1:10, kappa)
lines(1:10, kappa.est[2,])
lines(1:10, kappa.est[1,], lty = 2, col = "blue")
lines(1:10, kappa.est[3,], lty = 2, col = "blue")

plot.ts(kappa.chain[bn:it ,1])
plot.ts(kappa.chain[bn:it ,4])
plot.ts(kappa.chain[bn:it ,7])
plot.ts(kappa.chain[bn:it ,10])

acf(kappa.chain[bn:it ,1])
acf(kappa.chain[bn:it ,4])
acf(kappa.chain[bn:it ,7])
acf(kappa.chain[bn:it ,10])
graphics.off()


plot.ts(sigma_t.chain[bn:it])
abline(h = sigma_t, col = 2)

par(mfrow = c(1,3))
plot.ts(theta.chain[bn:it, 1])
abline(h = theta[1], col = 2)
plot.ts(theta.chain[bn:it, 2])
abline(h = theta[2], col = 2)
plot.ts(theta.chain[bn:it, 3])
abline(h = theta[3], col = 2)

### estimation (n-steps) #age-time-spatial parameters (no lambda, no  sp error) (after filtering) ----
it = 10000 #iterations
kappa.chain <- matrix(NA, nrow = it, ncol = t)

alpha.chain <- matrix(NA, nrow = it, ncol = ages)
beta.chain <- matrix(NA, nrow = it, ncol = ages)
gamma.chain <- matrix(NA, nrow = it, ncol = ages)

theta.chain <- matrix(NA, nrow = it, ncol = n)

sigma_w.chain <- rep(NA, it)
sigma_e.chain <- rep(NA, it)
#sigma_t.chain <- rep(NA, it)
#lambda.chain <- rep(NA, it)
nu.chain <- rep(NA, it)

## starters
alpha.chain[1, ] <- runif(ages)
beta.chain[1, ] <- runif(ages)
beta.chain[1, ] <- beta.chain[1, ]/sum(beta.chain[1, ])
gamma.chain[1, ] <- runif(ages)

theta.chain[1, ] <- theta    ### theta = 1 1 1 causou instabilidade

sigma_w.chain[1] <- 1
sigma_e.chain[1] <- 1
#sigma_t.chain[1] <- 1

#labda.min <- min(labda.dom$values); labda.max <- max(labda.dom$values)
#lambda.chain[1] <- runif(1, 1/labda.min, 1/labda.max)

nu.chain[1] <- runif(1, -1, 1)

## gibbs
start <- proc.time()[[3]]
set.seed(16)
for(i in 2:it){
  
  ### kappa
  fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                 sigma_e = sigma_e.chain[i-1],
                 sigma_w = sigma_w.chain[i-1],
                 nu = nu.chain[i-1],
                 alpha = alpha.chain[i-1, ],
                 beta = beta.chain[i-1, ],
                 gamma = gamma.chain[i-1, ],
                 theta = theta.chain[i-1, ])
  #theta = theta)
  
  kappa.chain[i,] <- fit$kappa.bs
  
  ### kappa constraint
  kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
  
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
    mean.reg <- aux.reg %*% t(X) %*% Y.aux
    var.reg <- sigma_e.chain[i] * aux.reg
    
    tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
    
    alpha.chain[i, j] <- tmp[1]
    beta.chain[i, j] <- tmp[2]
    gamma.chain[i, j] <- tmp[3]
  }
  
  ### alpha, beta & gamma constraint
  #alpha.chain[i, ] <- alpha.chain[i, ]/sum(alpha.chain[i, ])
  beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])
  #gamma.chain[i, ] <- gamma.chain[i, ]/sum(gamma.chain[i, ])
  
  ### sigma_w
  aux <- sum( (kappa.chain[i, 2:t] - nu.chain[i-1] - kappa.chain[i, 1:(t-1)])^2 ) ## years
  sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### nu
  aux <- sigma_w.chain[i]/t            #(t-1)
  aux.2 <- (kappa.chain[i, t] - kappa.chain[i, 1])/t  #(t-1)
  nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
  
  ### sigma_t
  # aux <- Mi - lambda*Wi
  # aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
  # #aux.2 <- theta %*% aux %*% matrix(theta)
  # sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
  ### um pouco esquisito valores mt altos surgirem da estimacao
  
  ### theta_q     ## rever contas, constraint nao deu jeito
  for(j in 1:n){
    mu0 = sum( (Wi[j, ]/Mi[j,j])*theta.chain[i-1, ] )
    sigma0 = sigma_t/Mi[j,j]
    
    y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta) %*% kappa.chain[i,]) - matrix(alpha,
                                                                                         nrow = ages,
                                                                                         ncol = t) )
    
    aux = ( sigma0*(t*sum(gamma^2)) + sigma_e )
    aux.2 = y.hat*matrix(gamma, nrow = ages, ncol = t)
    mu1 = ( sigma0*sum(aux.2) + sigma_e*mu0 )/aux
    sigma1 = (sigma_e*sigma0)/aux
    
    theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
  }
  
  ### theta constraint
  #theta.chain[i,] <- theta.chain[i,] - mean(theta.chain[i,])
  
  ### lambda
  
  
}
end <- proc.time()[[3]]
end - start ## seconds

### chain treatment
bn = 5000
thin = 5
kappa.est <- apply(kappa.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

alpha.est <- apply(alpha.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
beta.est <- apply(beta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
gamma.est <- apply(gamma.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
theta.est <- apply(theta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

sigma_w.est <- quantile(sigma_w.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
#sigma_t.est <- quantile(sigma_t.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_e.est <- quantile(sigma_e.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
nu.est <- quantile(nu.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))

#### age-time parameters plot
x11()
par(mfrow = c(2,2))
plot(1:10, kappa)
lines(1:10, kappa.est[2,])
lines(1:10, kappa.est[1,], lty = 2, col = "blue")
lines(1:10, kappa.est[3,], lty = 2, col = "blue")

plot(1:20, alpha)
lines(1:20, alpha.est[2,])
lines(1:20, alpha.est[1,], lty = 2, col = "blue")
lines(1:20, alpha.est[3,], lty = 2, col = "blue")

plot(1:20, beta)
lines(1:20, beta.est[2,])
lines(1:20, beta.est[1,], lty = 2, col = "blue")
lines(1:20, beta.est[3,], lty = 2, col = "blue")

plot(1:20, gamma)
lines(1:20, gamma.est[2,])
lines(1:20, gamma.est[1,], lty = 2, col = "blue")
lines(1:20, gamma.est[3,], lty = 2, col = "blue")

plot.ts(kappa.chain[bn:it ,1])
plot.ts(kappa.chain[bn:it ,4])
plot.ts(kappa.chain[bn:it ,7])
plot.ts(kappa.chain[bn:it ,10])

acf(kappa.chain[bn:it ,1])
acf(kappa.chain[bn:it ,4])
acf(kappa.chain[bn:it ,7])
acf(kappa.chain[bn:it ,10])

plot.ts(alpha.chain[bn:it ,1])
plot.ts(alpha.chain[bn:it ,8])
plot.ts(alpha.chain[bn:it ,15])
plot.ts(alpha.chain[bn:it ,20])

acf(alpha.chain[bn:it ,1])
acf(alpha.chain[bn:it ,8])
acf(alpha.chain[bn:it ,15])
acf(alpha.chain[bn:it ,20])

plot.ts(beta.chain[bn:it, 1])
plot.ts(beta.chain[bn:it, 7])
plot.ts(beta.chain[bn:it, 16])
plot.ts(beta.chain[bn:it, 20])

acf(beta.chain[bn:it, 1])
acf(beta.chain[bn:it, 7])
acf(beta.chain[bn:it, 16])
acf(beta.chain[bn:it, 20])

plot.ts(gamma.chain[bn:it, 1])
plot.ts(gamma.chain[bn:it, 4])
plot.ts(gamma.chain[bn:it, 12])
plot.ts(gamma.chain[bn:it, 20])

acf(gamma.chain[bn:it, 1])
acf(gamma.chain[bn:it, 4])
acf(gamma.chain[bn:it, 12])
acf(gamma.chain[bn:it, 20])
graphics.off()


nu.est
nu
sigma_w.est
sigma_w
sigma_e.est
sigma_e
theta.est
theta

par(mfrow = c(1, 3))
plot.ts(nu.chain[seq(bn, it, by = thin)])
abline(h = nu, col = 2)
plot.ts(sigma_w.chain[seq(bn, it, by = thin)])
abline(h = sigma_w, col = 2)
plot.ts(sigma_e.chain[seq(bn, it, by = thin)])
abline(h = sigma_e, col = 2)
#plot.ts(sigma_t.chain[bn:it])
#abline(h = sigma_t, col = 2)

par(mfrow = c(1,3))
plot.ts(theta.chain[seq(bn, it, by = thin), 1])
abline(h = theta[1], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 2])
abline(h = theta[2], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 3])
abline(h = theta[3], col = 2)

acf(nu.chain[bn:it])
acf(sigma_w.chain[bn:it])
acf(sigma_e.chain[bn:it])

### estimation (n-steps) #age-time-spatial parameters (no lambda) (after filtering) ----
it = 10000 #iterations
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
beta.chain[1, ] <- runif(ages)
beta.chain[1, ] <- beta.chain[1, ]/sum(beta.chain[1, ])
gamma.chain[1, ] <- runif(ages)
#gamma.chain[1, ] <- gamma.chain[1, ]/sum(gamma.chain[1, ])

theta.chain[1, ] <- theta    ### theta = 1 1 1 causou instabilidade

sigma_w.chain[1] <- 1
sigma_e.chain[1] <- 1
sigma_t.chain[1] <- 1

#labda.min <- min(labda.dom$values); labda.max <- max(labda.dom$values)
#lambda.chain[1] <- runif(1, 1/labda.min, 1/labda.max)

nu.chain[1] <- runif(1, -1, 1)

## gibbs
start <- proc.time()[[3]]
set.seed(16)
for(i in 2:it){
  
  ### kappa
  fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                 sigma_e = sigma_e.chain[i-1],
                 sigma_w = sigma_w.chain[i-1],
                 nu = nu.chain[i-1],
                 alpha = alpha.chain[i-1, ],
                 beta = beta.chain[i-1, ],
                 gamma = gamma.chain[i-1, ],
                 theta = theta.chain[i-1, ])
                 #theta = theta)
  
  kappa.chain[i,] <- fit$kappa.bs
  
  ### kappa constraint
  kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
  
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
    mean.reg <- aux.reg %*% t(X) %*% Y.aux
    var.reg <- sigma_e.chain[i] * aux.reg
    
    tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
    
    alpha.chain[i, j] <- tmp[1]
    beta.chain[i, j] <- tmp[2]
    gamma.chain[i, j] <- tmp[3]
  }
  
  ### alpha, beta & gamma constraint
  #alpha.chain[i, ] <- alpha.chain[i, ]/sum(alpha.chain[i, ])
  beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])
  #gamma.chain[i, ] <- gamma.chain[i, ]/sum(gamma.chain[i, ])
  
  ### sigma_w
  aux <- sum( (kappa.chain[i, 2:t] - nu.chain[i-1] - kappa.chain[i, 1:(t-1)])^2 ) ## years
  sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
  ### informative priors - gamma(a, b), a = 0.01, b = 0.01
  
  ### nu
  aux <- sigma_w.chain[i]/t            #(t-1)
  aux.2 <- (kappa.chain[i, t] - kappa.chain[i, 1])/t  #(t-1)
  nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
  
  ### sigma_t
  aux <- Mi - lambda*Wi
  aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
  #aux.2 <- theta %*% aux %*% matrix(theta)
  sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
  ### um pouco esquisito valores mt altos surgirem da estimacao
  
  ### theta_q     ## rever contas, constraint nao deu jeito
  for(j in 1:n){
    mu0 = sum( (Wi[j, ]/Mi[j,j])*theta.chain[i-1, ] )
    sigma0 = sigma_t/Mi[j,j]
    
    y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta) %*% kappa.chain[i,]) - matrix(alpha,
                                                                                         nrow = ages,
                                                                                         ncol = t) )
    
    aux = ( sigma0*(t*sum(gamma^2)) + sigma_e )
    aux.2 = y.hat*matrix(gamma, nrow = ages, ncol = t)
    mu1 = ( sigma0*sum(aux.2) + sigma_e*mu0 )/aux
    sigma1 = (sigma_e*sigma0)/aux
    
    theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
  }
  
  ### theta constraint
  #theta.chain[i,] <- theta.chain[i,] - mean(theta.chain[i,])
  
  # mu0 = lambda * rowSums(matrix(Wi/diag(Mi), nrow = n) %*% diag(theta.chain[i-1, ]))
  # sigma0 = sigma_t/diag(Mi)
  # y.hat = (y_obs - (beta.chain[i,] %*% t(kappa.chain[i,]) + alpha.chain[i,]))
  ### lambda
  
  
}
end <- proc.time()[[3]]
end - start ## seconds

### chain treatment
bn = 5000
thin = 2
kappa.est <- apply(kappa.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

alpha.est <- apply(alpha.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
beta.est <- apply(beta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
gamma.est <- apply(gamma.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
theta.est <- apply(theta.chain[seq(bn, it, by = thin), ], 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)

sigma_w.est <- quantile(sigma_w.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_t.est <- quantile(sigma_t.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
sigma_e.est <- quantile(sigma_e.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))
nu.est <- quantile(nu.chain[seq(bn, it, by = thin)], c(0.025, 0.5, 0.975))

#### age-time parameters plot
x11()
par(mfrow = c(2,2))
plot(1:10, kappa)
lines(1:10, kappa.est[2,])
lines(1:10, kappa.est[1,], lty = 2, col = "blue")
lines(1:10, kappa.est[3,], lty = 2, col = "blue")

plot(1:20, alpha)
lines(1:20, alpha.est[2,])
lines(1:20, alpha.est[1,], lty = 2, col = "blue")
lines(1:20, alpha.est[3,], lty = 2, col = "blue")

plot(1:20, beta)
lines(1:20, beta.est[2,])
lines(1:20, beta.est[1,], lty = 2, col = "blue")
lines(1:20, beta.est[3,], lty = 2, col = "blue")

plot(1:20, gamma)
lines(1:20, gamma.est[2,])
lines(1:20, gamma.est[1,], lty = 2, col = "blue")
lines(1:20, gamma.est[3,], lty = 2, col = "blue")

plot.ts(kappa.chain[bn:it ,1])
plot.ts(kappa.chain[bn:it ,4])
plot.ts(kappa.chain[bn:it ,7])
plot.ts(kappa.chain[bn:it ,10])

acf(kappa.chain[bn:it ,1])
acf(kappa.chain[bn:it ,4])
acf(kappa.chain[bn:it ,7])
acf(kappa.chain[bn:it ,10])

plot.ts(alpha.chain[bn:it ,1])
plot.ts(alpha.chain[bn:it ,8])
plot.ts(alpha.chain[bn:it ,15])
plot.ts(alpha.chain[bn:it ,20])

acf(alpha.chain[bn:it ,1])
acf(alpha.chain[bn:it ,8])
acf(alpha.chain[bn:it ,15])
acf(alpha.chain[bn:it ,20])

plot.ts(beta.chain[bn:it, 1])
plot.ts(beta.chain[bn:it, 7])
plot.ts(beta.chain[bn:it, 16])
plot.ts(beta.chain[bn:it, 20])

acf(beta.chain[bn:it, 1])
acf(beta.chain[bn:it, 7])
acf(beta.chain[bn:it, 16])
acf(beta.chain[bn:it, 20])

plot.ts(gamma.chain[bn:it, 1])
plot.ts(gamma.chain[bn:it, 4])
plot.ts(gamma.chain[bn:it, 12])
plot.ts(gamma.chain[bn:it, 20])

acf(gamma.chain[bn:it, 1])
acf(gamma.chain[bn:it, 4])
acf(gamma.chain[bn:it, 12])
acf(gamma.chain[bn:it, 20])
graphics.off()


nu.est
sigma_w.est
sigma_e.est
sigma_t.est
theta.est

#par(mfrow = c(1, 3))
plot.ts(nu.chain[seq(bn, it, by = thin)])
abline(h = nu, col = 2)
plot.ts(sigma_w.chain[seq(bn, it, by = thin)])
abline(h = sigma_w, col = 2)
plot.ts(sigma_e.chain[seq(bn, it, by = thin)])
abline(h = sigma_e, col = 2)
plot.ts(sigma_t.chain[seq(bn, it, by = thin)], ylim = c(0, 10))
abline(h = sigma_t, col = 2)

par(mfrow = c(1,3))
plot.ts(theta.chain[seq(bn, it, by = thin), 1])
abline(h = theta[1], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 2])
abline(h = theta[2], col = 2)
plot.ts(theta.chain[seq(bn, it, by = thin), 3])
abline(h = theta[3], col = 2)

acf(nu.chain[bn:it])
acf(sigma_w.chain[bn:it])
acf(sigma_e.chain[bn:it])
acf(sigma_t.chain[bn:it])
graphics.off()
