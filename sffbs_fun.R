
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

