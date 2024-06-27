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
Q = M - 1*W

### assumindo apenas 3 microreg pra simplificar:
ages <- 20  ## idades
t <- 10  ## tempo
n <- 3   ## regioes

Qi <- Q[1:n, 1:n]
set.seed(14)
theta <- MASS::mvrnorm(n = 1, mu = rep(0, n), Sigma = MASS::ginv(Qi))


### alpha, beta and kappa ----
set.seed(14)
alpha <- runif(ages, -1, 1)
beta <- runif(ages, 0, 1)
kappa0 <- 1
kappa <- rep(NA, t)
nu <- 1
set.seed(14)
kappa[1] <- nu + kappa0 + rnorm(1, 0, sqrt(0.025)) 
for(k in 2:t){
  kappa[k] <- nu + kappa[k-1] + rnorm(1, 0, sqrt(0.025))
}
gamma <- rep(1, ages)

### aux ----
A <- kronecker(rep(1, n), alpha)
B <- kronecker(rep(1, n), beta)
C <- kronecker(theta, gamma)
set.seed(14)
y_obs <- matrix(NA, nrow = ages*n, ncol = t)
for(k in 1:t){
  y_obs[,k] <- A + (B * kappa[k]) + C + rnorm(n*ages, 0, sd = sqrt(0.05))
}


### estimation ----
#### FF
m0 <- 0; C0 <- 100
sigma_e = 0.05
sigma_w = 0.025
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


plot(1:10, kappa)
lines(1:10, mt)

kappa.ff <- rnorm(60, mean = at, sd = sqrt(Rt))

par(mfrow=c(1,3))
plot(1:20, mt[1:20], type = "l", main = "m = 1")
points(1:20, y_obs[1:20])
plot(21:40, mt[21:40], type = "l", main = "m = 2")
points(21:40, y_obs[21:40])
plot(41:60, mt[41:60], type = "l", main = "m = 3")
points(41:60, y_obs[41:60])

kappa.bs <- rep(NA, t)
set.seed(14)
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

plot(1:10, kappa)
lines(1:10, kappa.bs)

### medias (?)
par(mfrow=c(2,3))
plot(1:20, mt[1:20], type = "l", main = "m.ff = 1")
points(1:20, y_obs[1:20])
plot(21:40, mt[21:40], type = "l", main = "m.ff = 2")
points(21:40, y_obs[21:40])
plot(41:60, mt[41:60], type = "l", main = "m.ff = 3")
points(41:60, y_obs[41:60])

plot(1:20, s0[1:20], type = "l", main = "m.bs = 1")
points(1:20, y_obs[1:20])
plot(21:40, s0[21:40], type = "l", main = "m.bs = 2")
points(21:40, y_obs[21:40])
plot(41:60, s0[41:60], type = "l", main = "m.bs = 3")
points(41:60, y_obs[41:60])

### kappas (?)
par(mfrow=c(2,3))
plot(1:20, kappa.ff[1:20], type = "l", main = "k.ff = 1")
points(1:20, kappa[1:20])
plot(21:40, kappa.ff[21:40], type = "l", main = "k.ff = 2")
points(21:40, kappa[21:40])
plot(41:60, kappa.ff[41:60], type = "l", main = "k.ff = 3")
points(41:60, kappa[41:60])

plot(1:20, kappa.bs[1:20], type = "l", main = "k.bs = 1")
points(1:20, kappa[1:20])
plot(21:40, kappa.bs[21:40], type = "l", main = "k.bs = 2")
points(21:40, kappa[21:40])
plot(41:60, kappa.bs[41:60], type = "l", main = "k.bs = 3")
points(41:60, kappa[41:60])

