
sblc <- function(y, ages, t, n, m0, C0, M, W, it, bn, thin){
  
  y_obs = y
  
  ### CHAINS
  kappa.chain <- matrix(NA, nrow = it, ncol = t+1)
  
  alpha.chain <- matrix(NA, nrow = it, ncol = ages)
  beta.chain <- matrix(NA, nrow = it, ncol = ages)
  gamma.chain <- matrix(NA, nrow = it, ncol = ages)
  
  theta.chain <- matrix(NA, nrow = it, ncol = n)
  
  sigma_w.chain <- rep(NA, it)
  sigma_e.chain <- rep(NA, it)
  sigma_t.chain <- rep(NA, it)
  nu.chain <- rep(NA, it)
  
  ### STARTERS
  alpha.chain[1, ] <- runif(ages)
  beta.chain[1, ] <- runif(ages)
  beta.chain[1, ] <- beta.chain[1, ]/sum(beta.chain[1, ])
  gamma.chain[1, ] <- runif(ages)
  gamma.chain[1, ] <- gamma.chain[1, ]/sum(gamma.chain[1, ])
  
  theta.chain[1, ] <- runif(n)
  theta.chain[1, ] <- theta.chain[1, ] - mean(theta.chain[1, ])
  
  sigma_w.chain[1] <- 1
  sigma_e.chain[1] <- 1
  sigma_t.chain[1] <- 1
  nu.chain[1] <- runif(1, -1, 1)
  
  ### GIBBS
  pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = it, clear = FALSE, width = 60)
  pb$tick()
  prop = 0
  
  set.seed(16)
  for(i in 2:it){
    
    pb$tick()
    
    ### KAPPA
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
    ### KAPPA CONSTRAINT
    kappa.chain[i,] <- kappa.chain[i, ] - mean(kappa.chain[i, ])
    
    ### SIGMA_E
    A <- matrix(kronecker(rep(1, n), alpha.chain[i-1, ]))
    B <- matrix(kronecker(rep(1, n), beta.chain[i-1, ]))
    C <- matrix(kronecker(theta.chain[i-1, ], gamma.chain[i-1, ]))
    
    aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i, 2:(t+1)]) - C[, rep(1, t)])^2 )
    sigma_e.chain[i] <- invgamma::rinvgamma(1, t*n*ages/2, rate = aux/2) ##gerar fixando idade?
    ### informative priors - gamma(a, b), a = 0.01, b = 0.01
    
    ### ALPHA, BETA AND GAMMA
    for (j in 1:ages) {
      Y.aux <- c( y_obs[j, ] )
      for(k in 2:n){
        Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
      }
      
      X <- cbind(1, rep(kappa.chain[i, 2:(t+1)], n), rep(theta.chain[i-1, ], each = t))
      aux.reg <- chol2inv(chol(t(X) %*% X))
      mean.reg <- aux.reg %*% t(X) %*% matrix(Y.aux)
      var.reg <- sigma_e.chain[i] * aux.reg
      
      tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
      
      alpha.chain[i, j] <- tmp[1]
      beta.chain[i, j] <- tmp[2]
      gamma.chain[i, j] <- tmp[3]
    }
    
    ### BETA & GAMMA CONSTRAINT
    beta.chain[i, ] <- beta.chain[i, ]/sum(beta.chain[i, ])
    gamma.chain[i, ] <- gamma.chain[i, ]/sum(gamma.chain[i, ])
    
    ### SIGMA_W
    aux <- sum( (kappa.chain[i, 2:(t+1)] - nu.chain[i-1] - kappa.chain[i, 1:t])^2 ) ## years
    sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
    ### informative priors - gamma(a, b), a = 0.01, b = 0.01
    
    ### NU
    aux <- sigma_w.chain[i]/(t-1)
    aux.2 <- (kappa.chain[i, t+1] - kappa.chain[i, 1])/t
    nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
    #nu.chain[i] <- nu
    
    ### SIGMA_T
    aux <- M - W
    aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
    sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
    
    ### THETA_q
    for(j in 1:n){
      mu0 = sum( (W[j, -j]/M[j,j])*theta.chain[i-1, -j] )
      sigma0 = sigma_t.chain[i]/M[j,j]

      
      y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta.chain[i,]) %*% kappa.chain[i, 2:(t+1)]) - matrix(alpha.chain[i,],
                                                                                                             nrow = ages,
                                                                                                             ncol = t) )
      
      aux = ( sigma0*(t*sum(gamma.chain[i,]^2)) + sigma_e.chain[i] )
      aux.2 = y.hat*matrix(gamma.chain[i,], nrow = ages, ncol = t)
      mu1 = ( sigma0*sum(aux.2) + sigma_e.chain[i]*mu0 )/aux
      sigma1 = (sigma_e.chain[i]*sigma0)/aux
      
      theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
    }
    
    ### THETA CONSTRAINT
    theta.chain[i,] <- theta.chain[i,] - mean(theta.chain[i,])

  }
  
  ### CHAIN TREATMENT
  kappa.est <- kappa.chain[seq(bn, it, by = thin), 2:(t+1)]
  kappa0.est <- kappa.chain[seq(bn, it, by = thin), 1]
  
  alpha.est <- alpha.chain[seq(bn, it, by = thin), ]
  beta.est <- beta.chain[seq(bn, it, by = thin), ]
  gamma.est <- gamma.chain[seq(bn, it, by = thin), ]
  theta.est <- theta.chain[seq(bn, it, by = thin), ]
  
  sigma_w.est <- sigma_w.chain[seq(bn, it, by = thin)]
  sigma_t.est <- sigma_t.chain[seq(bn, it, by = thin)]
  sigma_e.est <- sigma_e.chain[seq(bn, it, by = thin)]
  nu.est <- nu.chain[seq(bn, it, by = thin)]
  
  
  ### RETURN
  fit <- list(info = list(y = y_obs,
                          ages = ages,
                          t = t,
                          n = n,
                          m0 = m0,
                          C0 = C0,
                          M = M,
                          W = W,
                          it = it,
                          bn = bn,
                          thin = thin),
              kappa.chain = kappa.est,
              kappa0.chain = kappa0.est,
              alpha.chain = alpha.est,
              beta.chain = beta.est,
              gamma.chain = gamma.est,
              theta.chain = theta.est,
              sigma_w.chain = sigma_w.est,
              sigma_e.chain = sigma_e.est,
              sigma_t.chain = sigma_t.est,
              nu.chain = nu.est)
}


sblc_2 <- function(y, ages, t, n, m0, C0, M, W, it, bn, thin){
  
  y_obs = y
  
  ### CHAINS
  #kappa.chain <- matrix(NA, nrow = it, ncol = t+1)
  kappa.chain <- matrix(NA, nrow = it, ncol = t)
  
  alpha.chain <- matrix(NA, nrow = it, ncol = ages)
  beta.chain <- matrix(NA, nrow = it, ncol = ages)
  gamma.chain <- matrix(NA, nrow = it, ncol = ages)
  
  theta.chain <- matrix(NA, nrow = it, ncol = n)
  
  sigma_w.chain <- rep(NA, it)
  sigma_e.chain <- rep(NA, it)
  sigma_t.chain <- rep(NA, it)
  nu.chain <- rep(NA, it)
  
  ### STARTERS
  alpha.chain[1, ] <- runif(ages)
  beta.chain[1, ] <- runif(ages)
  #beta.chain[1, ] <- beta.chain[1, ]/sum(beta.chain[1, ])
  gamma.chain[1, ] <- runif(ages)
  #gamma.chain[1, ] <- gamma.chain[1, ]/sum(gamma.chain[1, ])
  
  theta.chain[1, ] <- runif(n)
  #theta.chain[1, ] <- theta.chain[1, ] - mean(theta.chain[1, ])
  
  sigma_w.chain[1] <- 1
  sigma_e.chain[1] <- 1
  sigma_t.chain[1] <- 1
  nu.chain[1] <- runif(1, -1, 1)
  
  ### GIBBS
  pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = it, clear = FALSE, width = 60)
  pb$tick()
  prop = 0
  
  set.seed(16)
  for(i in 2:it){
    
    pb$tick()
    
    ### KAPPA
    fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                   sigma_e = sigma_e.chain[i-1],
                   sigma_w = sigma_w.chain[i-1],
                   nu = nu.chain[i-1],
                   alpha = alpha.chain[i-1, ],
                   beta = beta.chain[i-1, ],
                   gamma = gamma.chain[i-1, ],
                   theta = theta.chain[i-1, ])
    
    # kappa.chain[i, 2:(t+1)] <- fit$kappa.bs
    # kappa.chain[i, 1] <- fit$kappa0.bs
    kappa.chain[i, ] <- fit$kappa.bs
    
    
    ### CONSTRAINTS
    level <- mean(kappa.chain[i,])
    incline <- 1/sum(beta.chain[i-1,])
    #incline <- (chain$kappa[1,i] - chain$kappa[N,i]) / (N - 1)
    
    kappa.chain[i,] <- (kappa.chain[i, ] - level) / incline
    
    level2 <- mean(theta.chain[i-1,])
    incline2 <- 1/sum(gamma.chain[i-1,])
    
    theta.chain[i-1,] <- (theta.chain[i-1,] - level2) / incline2
    
    alpha.chain[i-1,] <- alpha.chain[i-1,] + (beta.chain[i-1,] * level) + (gamma.chain[i-1,] * level2)
    beta.chain[i-1,] <- beta.chain[i-1,] * incline
    gamma.chain[i-1,] <- gamma.chain[i-1,] * incline2
    
  
    
    ### SIGMA_E
    A <- matrix(kronecker(rep(1, n), alpha.chain[i-1, ]))
    B <- matrix(kronecker(rep(1, n), beta.chain[i-1, ]))
    C <- matrix(kronecker(theta.chain[i-1, ], gamma.chain[i-1, ]))
    
    aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i, ]) - C[, rep(1, t)])^2 )
    sigma_e.chain[i] <- invgamma::rinvgamma(1, t*n*ages/2, rate = aux/2) ##gerar fixando idade?
    ### informative priors - gamma(a, b), a = 0.01, b = 0.01
    
    ### ALPHA, BETA AND GAMMA
    for (j in 1:ages) {
      Y.aux <- c( y_obs[j, ] )
      for(k in 2:n){
        Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
      }
      
      X <- cbind(1, rep(kappa.chain[i, ], n), rep(theta.chain[i-1, ], each = t))
      aux.reg <- chol2inv(chol(t(X) %*% X))
      mean.reg <- aux.reg %*% t(X) %*% matrix(Y.aux)
      var.reg <- sigma_e.chain[i] * aux.reg
      
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
    aux <- M - W
    aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
    sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
    
    ### THETA_q
    for(j in 1:n){
      mu0 = sum( (W[j, -j]/M[j,j])*theta.chain[i-1, -j] )
      sigma0 = sigma_t.chain[i]/M[j,j]
      
      
      y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta.chain[i,]) %*% kappa.chain[i, ]) - matrix(alpha.chain[i,],
                                                                                                             nrow = ages,
                                                                                                             ncol = t) )
      
      aux = ( sigma0*(t*sum(gamma.chain[i,]^2)) + sigma_e.chain[i] )
      aux.2 = y.hat*matrix(gamma.chain[i,], nrow = ages, ncol = t)
      mu1 = ( sigma0*sum(aux.2) + sigma_e.chain[i]*mu0 )/aux
      sigma1 = (sigma_e.chain[i]*sigma0)/aux
      
      theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
    }
    
    ### THETA CONSTRAINT
    # theta.chain[i,] <- theta.chain[i,] - mean(theta.chain[i,])
    
  }
  
  ### CHAIN TREATMENT
  kappa.est <- kappa.chain[seq(bn, it, by = thin), ]
  #kappa0.est <- kappa.chain[seq(bn, it, by = thin), 1]
  
  alpha.est <- alpha.chain[seq(bn, it, by = thin), ]
  beta.est <- beta.chain[seq(bn, it, by = thin), ]
  gamma.est <- gamma.chain[seq(bn, it, by = thin), ]
  theta.est <- theta.chain[seq(bn, it, by = thin), ]
  
  sigma_w.est <- sigma_w.chain[seq(bn, it, by = thin)]
  sigma_t.est <- sigma_t.chain[seq(bn, it, by = thin)]
  sigma_e.est <- sigma_e.chain[seq(bn, it, by = thin)]
  nu.est <- nu.chain[seq(bn, it, by = thin)]
  
  
  ### RETURN
  fit <- list(info = list(y = y_obs,
                          ages = ages,
                          t = t,
                          n = n,
                          m0 = m0,
                          C0 = C0,
                          M = M,
                          W = W,
                          it = it,
                          bn = bn,
                          thin = thin),
              kappa.chain = kappa.est,
              #kappa0.chain = kappa0.est,
              alpha.chain = alpha.est,
              beta.chain = beta.est,
              gamma.chain = gamma.est,
              theta.chain = theta.est,
              sigma_w.chain = sigma_w.est,
              sigma_e.chain = sigma_e.est,
              sigma_t.chain = sigma_t.est,
              nu.chain = nu.est)
}


sblc_missing <- function(y, ages, t, n, m0, C0, M, W, it, bn, thin){
  
  y_obs = y
  ind_missing = which(is.na(y), arr.ind = T) ## missing index
  
  ### CHAINS
  #kappa.chain <- matrix(NA, nrow = it, ncol = t+1)
  kappa.chain <- matrix(NA, nrow = it, ncol = t)
  
  alpha.chain <- matrix(NA, nrow = it, ncol = ages)
  beta.chain <- matrix(NA, nrow = it, ncol = ages)
  gamma.chain <- matrix(NA, nrow = it, ncol = ages)
  
  theta.chain <- matrix(NA, nrow = it, ncol = n)
  
  sigma_w.chain <- rep(NA, it)
  sigma_e.chain <- rep(NA, it)
  sigma_t.chain <- rep(NA, it)
  nu.chain <- rep(NA, it)
  
  input.chain <- matrix(NA, nrow = it, ncol = nrow(ind_missing))
  
  ### STARTERS
  alpha.chain[1, ] <- runif(ages)
  beta.chain[1, ] <- runif(ages)
  gamma.chain[1, ] <- runif(ages)
  
  theta.chain[1, ] <- runif(n)
  
  sigma_w.chain[1] <- 1
  sigma_e.chain[1] <- 1
  sigma_t.chain[1] <- 1
  nu.chain[1] <- runif(1, -1, 1)
  
  input.chain[1,] <- rowMeans(y, na.rm = T)[ind_missing[,1]]
  y_obs[ind_missing] <- input.chain[1,]
  
  ### GIBBS
  pb  = progress::progress_bar$new(format = "Simulating [:bar] :percent in :elapsed",total = it, clear = FALSE, width = 60)
  pb$tick()
  prop = 0
  
  set.seed(16)
  for(i in 2:it){
    
    pb$tick()
    
    ### KAPPA
    fit <- ffbs_sp(y_obs = y_obs, ages = ages, t = t, n = n, m0 = m0, C0 = C0,
                   sigma_e = sigma_e.chain[i-1],
                   sigma_w = sigma_w.chain[i-1],
                   nu = nu.chain[i-1],
                   alpha = alpha.chain[i-1, ],
                   beta = beta.chain[i-1, ],
                   gamma = gamma.chain[i-1, ],
                   theta = theta.chain[i-1, ])
    
    # kappa.chain[i, 2:(t+1)] <- fit$kappa.bs
    # kappa.chain[i, 1] <- fit$kappa0.bs
    kappa.chain[i, ] <- fit$kappa.bs
    
    
    ### CONSTRAINTS
    level <- mean(kappa.chain[i,])
    incline <- sum(beta.chain[i-1,])
    #incline <- (kappa.chain[i,1] - kappa.chain[i,t]) / (t - 1)
    level2 <- mean(theta.chain[i-1,])
    incline2 <- sum(gamma.chain[i-1,])
    
    
    alpha.chain[i-1,] <- alpha.chain[i-1,] + (beta.chain[i-1,] * level) + (gamma.chain[i-1,] * level2)
    kappa.chain[i,] <- (kappa.chain[i, ] - level) * incline
    theta.chain[i-1,] <- (theta.chain[i-1,] - level2) * incline2
    beta.chain[i-1,] <- beta.chain[i-1,] / incline
    gamma.chain[i-1,] <- gamma.chain[i-1,] / incline2
    
    
    ### SIGMA_E
    A <- matrix(kronecker(rep(1, n), alpha.chain[i-1, ]))
    B <- matrix(kronecker(rep(1, n), beta.chain[i-1, ]))
    C <- matrix(kronecker(theta.chain[i-1, ], gamma.chain[i-1, ]))
    
    aux <- sum( (y_obs - A[, rep(1, t)] - (B %*% kappa.chain[i, ]) - C[, rep(1, t)])^2 )
    sigma_e.chain[i] <- invgamma::rinvgamma(1, t*n*ages/2, rate = aux/2) ##gerar fixando idade?
    ### informative priors - gamma(a, b), a = 0.01, b = 0.01
    
    ### ALPHA, BETA AND GAMMA
    for (j in 1:ages) {
      Y.aux <- c( y_obs[j, ] )
      for(k in 2:n){
        Y.aux <- c( Y.aux, y_obs[( j + (ages * (k-1)) ), ] )
      }
      
      X <- cbind(1, rep(kappa.chain[i, ], n), rep(theta.chain[i-1, ], each = t))
      aux.reg <- chol2inv(chol(t(X) %*% X))
      mean.reg <- aux.reg %*% t(X) %*% matrix(Y.aux)
      var.reg <- sigma_e.chain[i] * aux.reg
      
      tmp <- MASS::mvrnorm(1, mean.reg, var.reg)
      
      alpha.chain[i, j] <- tmp[1]
      beta.chain[i, j] <- tmp[2]
      gamma.chain[i, j] <- tmp[3]
    }
    
    ### SIGMA_W
    aux <- sum( (kappa.chain[i, 2:t] - nu.chain[i-1] - kappa.chain[i, 1:(t-1)])^2 ) ## years
    sigma_w.chain[i] <- invgamma::rinvgamma(1, (t-1)/2, rate = aux/2)
    ### informative priors - gamma(a, b), a = 0.01, b = 0.01
    
    ### NU
    aux <- sigma_w.chain[i]/(t-1)
    aux.2 <- (kappa.chain[i, t] - kappa.chain[i, 1])/(t-1)
    nu.chain[i] <- rnorm(1, aux.2, sqrt(aux))
    
    ### SIGMA_T
    aux <- M - W
    aux.2 <- theta.chain[i-1, ] %*% aux %*% matrix(theta.chain[i-1, ])
    sigma_t.chain[i] <- invgamma::rinvgamma(1, (n-1)/2, rate = aux.2/2)
    
    ### THETA_q
    for(j in 1:n){
      mu0 = sum( (W[j, -j]/M[j,j])*theta.chain[i-1, -j] )
      sigma0 = sigma_t.chain[i]/M[j,j]
      
      
      y.hat = ( y_obs[1:ages + ages*(j-1), ] - (matrix(beta.chain[i,]) %*% kappa.chain[i, ]) - matrix(alpha.chain[i,],
                                                                                                      nrow = ages,
                                                                                                      ncol = t) )
      
      aux = ( sigma0*(t*sum(gamma.chain[i,]^2)) + sigma_e.chain[i] )
      aux.2 = y.hat*matrix(gamma.chain[i,], nrow = ages, ncol = t)
      mu1 = ( sigma0*sum(aux.2) + sigma_e.chain[i]*mu0 )/aux
      sigma1 = (sigma_e.chain[i]*sigma0)/aux
      
      theta.chain[i,j] <- rnorm(1, mu1, sqrt(sigma1))
    }
    
    
    ### MISSING
    for(j in 1:nrow(ind_missing)){
      loc <- c(ind_missing[j,1] %% ages,
               ind_missing[j,2],
               ifelse(ind_missing[j,1]%%ages == 0, ind_missing[j,1]/ages, ind_missing[j,1]%/%ages + 1) )
      
      aux = alpha.chain[i, loc[1]] + (beta.chain[i, loc[1]] * kappa.chain[i, loc[2]]) + (gamma.chain[i, loc[1]] * theta.chain[i, loc[3]])
      input.chain[i,j] <- rnorm(1, mean = aux, sd = sqrt(sigma_e.chain[i]))
    }
    y_obs[ind_missing] <- input.chain[i,]
    
  }
  
  ### CHAIN TREATMENT
  kappa.est <- kappa.chain[seq(bn, it, by = thin), ]
  #kappa0.est <- kappa.chain[seq(bn, it, by = thin), 1]
  
  alpha.est <- alpha.chain[seq(bn, it, by = thin), ]
  beta.est <- beta.chain[seq(bn, it, by = thin), ]
  gamma.est <- gamma.chain[seq(bn, it, by = thin), ]
  theta.est <- theta.chain[seq(bn, it, by = thin), ]
  
  sigma_w.est <- sigma_w.chain[seq(bn, it, by = thin)]
  sigma_t.est <- sigma_t.chain[seq(bn, it, by = thin)]
  sigma_e.est <- sigma_e.chain[seq(bn, it, by = thin)]
  nu.est <- nu.chain[seq(bn, it, by = thin)]
  
  input.est <- input.chain[seq(bn, it, by = thin), ]
  
  
  ### RETURN
  fit <- list(info = list(y = y,
                          ages = ages,
                          t = t,
                          n = n,
                          m0 = m0,
                          C0 = C0,
                          M = M,
                          W = W,
                          it = it,
                          bn = bn,
                          thin = thin),
              kappa.chain = kappa.est,
              #kappa0.chain = kappa0.est,
              alpha.chain = alpha.est,
              beta.chain = beta.est,
              gamma.chain = gamma.est,
              theta.chain = theta.est,
              sigma_w.chain = sigma_w.est,
              sigma_e.chain = sigma_e.est,
              sigma_t.chain = sigma_t.est,
              nu.chain = nu.est,
              input.chain = input.est)
}
