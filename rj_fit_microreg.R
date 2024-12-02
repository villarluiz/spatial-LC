library(sf)
library(spdep)
library(tidyverse)

#rj_mun <- st_read("RJ_Municipios_2022/RJ_Municipios_2022.shp")
rj_micro <- st_read("RJ_Microrregioes_2022/RJ_Microrregioes_2022.shp")
#plot(rj_mun[,c(2,5)])
#plot(rj_micro[,c(2,5)])

shape_aisp_nb_val <- st_make_valid(rj_micro[,c(2,5)])
nb <- poly2nb(shape_aisp_nb_val, queen=TRUE)

# matriz de vizinhaça
W <- nb2mat(nb, style="B", zero.policy = TRUE)
M <- diag(rowSums(W))

ages <- 13  ## idades    20-24, 25-29, ..., 80+
t <- 42  ## tempo
n <- 18   ## regioes

it = 20000 #iterations
bn = 10000
thin = 4

m0 <- 0; C0 <- 100

### y_obs
y <- matrix(NA, nrow = ages*n, ncol = t)
for(i in 1:t){
  aux <- i + 1979
  mx <- as.data.frame(read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt")))[,-1] %>%
    stack()
  y[,i] <- mx$values
}
ind <- as.character(mx$ind)



### missing data input (interpolando no ano, nao no espaço)
y_treat <- t(zoo::na.approx(t(y)))
y_treat[92, 42] <- y_treat[92, 41] ## repetindo valor do ultimo ano (2021) pq nao tem como interpolar


### moran
result <- matrix(NA, nrow = ages, ncol = t)
data_w <- nb2listw(nb, style = "W")
for(k in 1:t){
  for(i in 1:ages){
    aux <- y_treat[seq(i, 234, by = 13), k]
    result[i,k] <- moran.test(aux, data_w)$statistic
  }
}

options(scipen=999)
library(tidyverse)
colnames(result) <- 1980:2021
result.2 <- cbind(as.data.frame(result), x = 1:13) %>%
  pivot_longer(1:42)
plotly::ggplotly(ggplot(result.2) +
                   scale_x_discrete(breaks = scales::pretty_breaks()) +
                   geom_tile(aes(x = name, y = x, fill = value)) +
                   scale_fill_gradient2(low = "red4", mid = "white", high = "blue", midpoint = 0) +
                   labs(x = "year", y = "age interval", fill = "p-value") )

sum(result > 0.05)/(ages*t)   ###78.5% dos anos/tempo tem correlacao 




source("sblc_fun.R")
source("sffbs_fun.R")
fit <- sblc(y_treat, ages, t, n, m0, C0, M, W, it = 2000, bn = 1000, thin = 1)

fit2 <- sblc_2(y_treat, ages, t, n, m0, C0, M, W, it = 2000, bn = 1000, thin = 1)
#saveRDS(fit2, "rj_micro_fit2.RDS")
#saveRDS(fit2, "rj_micro_fit2_full.RDS")

###ajuste sigma_e var   + var, + suav
source("sffbs_fun_sigmae.R")
fit3 <- sblc_3(y_treat, ages, t, n, m0, C0, M, W, it = 20000, bn = 10000, thin = 5)
#saveRDS(fit3, "rj_micro_fit3.RDS")
saveRDS(fit3, "rj_micro_fit3_full.RDS")

fit4 <- sblc_3_2(y_treat, ages, t, n, m0, C0, M, W, it = 2000, bn = 1000, thin = 1)

saveRDS(fit4, "rj_micro_fit3_2.RDS")

fit2 <- readRDS("rj_micro_fit2.RDS")
fit3 <- readRDS("rj_micro_fit3.RDS")
fit3.full <- readRDS("rj_micro_fit3_full.RDS")
fit = fit4
### vis. cadeias

plot.ts(fit$kappa.chain[ ,1])
plot.ts(fit$kappa.chain[ ,10])
plot.ts(fit$kappa.chain[ ,20])
plot.ts(fit$kappa.chain[ ,30])
plot.ts(fit$kappa.chain[ ,42])

plot.ts(fit2$kappa.chain[ ,1])
plot.ts(fit2$kappa.chain[ ,10])
plot.ts(fit2$kappa.chain[ ,20])
plot.ts(fit2$kappa.chain[ ,30])

## é alguma questão no ffbs, que os kappas estao com um drift mt pequeno, era pra ser -0.95

acf(fit$kappa.chain[ ,1])
acf(fit$kappa.chain[ ,10])
acf(fit$kappa.chain[ ,20])
acf(fit$kappa.chain[ ,30])

plot.ts(fit$alpha.chain[ ,1])
plot.ts(fit$alpha.chain[ ,10])
plot.ts(fit$alpha.chain[ ,13])

acf(fit$alpha.chain[ ,1])
acf(fit$alpha.chain[ ,10])
acf(fit$alpha.chain[ ,13])

plot.ts(fit$beta.chain[ ,1])
plot.ts(fit$beta.chain[ ,10])
plot.ts(fit$beta.chain[ ,13])

acf(fit$beta.chain[ ,1])
acf(fit$beta.chain[ ,10])
acf(fit$beta.chain[ ,13])

plot.ts(fit$gamma.chain[ ,1])
plot.ts(fit$gamma.chain[ ,10])
plot.ts(fit$gamma.chain[ ,13])

acf(fit$gamma.chain[ ,1])
acf(fit$gamma.chain[ ,10])
acf(fit$gamma.chain[ ,13])

plot.ts(fit$theta.chain[ ,1])
plot.ts(fit$theta.chain[ ,10])
plot.ts(fit$theta.chain[ ,13])
plot.ts(fit$theta.chain[ ,n])

acf(fit$theta.chain[ ,1])
acf(fit$theta.chain[ ,10])
acf(fit$theta.chain[ ,13])

graphics.off()

plot.ts(fit$nu.chain)
acf(fit$nu.chain)
plot.ts(fit$sigma_w.chain)  ### diminuiu ebaaa
acf(fit$sigma_w.chain)
plot.ts(fit$sigma_e.chain[,1])
plot.ts(fit$sigma_e.chain[,10])
plot.ts(fit$sigma_e.chain[,13])

plot.ts(fit$sigma_t.chain)  ##diminuiu
acf(fit$sigma_t.chain)

plot.ts(apply(fit$theta.chain, 2, median))
### ok sigma_t ser 5?


plot(apply(fit$kappa.chain, 2, median), type = "l")
lines(apply(fit$kappa.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$kappa.chain, 2, quantile, 0.975), col = 2, lty = 2)

lines(apply(fit3$kappa.chain, 2, median), col = 3)
lines(apply(fit3$kappa.chain, 2, quantile, 0.025), col = 3, lty = 2)
lines(apply(fit3$kappa.chain, 2, quantile, 0.975), col = 3, lty = 2)

plot(apply(fit$alpha.chain, 2, median), type = "l")
lines(apply(fit$alpha.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$alpha.chain, 2, quantile, 0.975), col = 2, lty = 2)
lines(apply(fit3$alpha.chain, 2, median), col = 3)
lines(apply(fit3$alpha.chain, 2, quantile, 0.025), col = 3, lty = 2)
lines(apply(fit3$alpha.chain, 2, quantile, 0.975), col = 3, lty = 2)

plot(apply(fit$beta.chain, 2, median), type = "l")
lines(apply(fit$beta.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$beta.chain, 2, quantile, 0.975), col = 2, lty = 2)
lines(apply(fit3$beta.chain, 2, median), col = 3)
lines(apply(fit3$beta.chain, 2, quantile, 0.025), col = 3, lty = 2)
lines(apply(fit3$beta.chain, 2, quantile, 0.975), col = 3, lty = 2)
abline(h = 0)

plot(apply(fit$gamma.chain, 2, median), type = "l")
lines(apply(fit$gamma.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$gamma.chain, 2, quantile, 0.975), col = 2, lty = 2)
lines(apply(fit3$gamma.chain, 2, median), col = 3)
lines(apply(fit3$gamma.chain, 2, quantile, 0.025), col = 3, lty = 2)
lines(apply(fit3$gamma.chain, 2, quantile, 0.975), col = 3, lty = 2)

plot(apply(fit$theta.chain, 2, median), type = "l")
lines(apply(fit$theta.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$theta.chain, 2, quantile, 0.975), col = 2, lty = 2)
lines(apply(fit3$theta.chain, 2, median), col = 3)
lines(apply(fit3$theta.chain, 2, quantile, 0.025), col = 3, lty = 2)
lines(apply(fit3$theta.chain, 2, quantile, 0.975), col = 3, lty = 2)
### maior intervalo (interessante)
### deslocado quando comparado com theta_n = 0

#### gerar y_ts
source("fitted_sblc.R")   ###funciona tambem pra sigmae variando na idade
fit <- readRDS("rj_micro_fit.RDS")

fit2.qx <- fitted_sblc(fit2)
fit3.qx <- fitted_sblc(fit3.full)
fit4.qx <- fitted_sblc(fit4)

## comparar qxs via metrica residual pearson quadrada = (o - e)^2/e, fazer grafico heatmap por reg
library(ggplot2)
### y_obs
ex <- matrix(NA, nrow = ages*n, ncol = t)
dx <- matrix(NA, nrow = ages*n, ncol = t)
for(i in 1:t){
  aux <- i + 1979
  mx <- as.data.frame(t(read.table(paste0("Dados_microrregioes/det/ex_det/ex_", as.character(aux), ".txt"), skip = 1, sep  =";")))[-c(1,15),-19] %>%
    stack()
  ex[,i] <- as.numeric(mx$values)
  mx <- as.data.frame(t(read.table(paste0("Dados_microrregioes/det/obt_det/obt_", as.character(aux), ".txt"), skip = 1, sep  =";")))[-c(1,15),-19] %>%
    stack()
  dx[,i] <- as.numeric(mx$values)
}


aux.qx.2 <- fit2.qx$mean
aux.qx.3 <- fit3.qx$mean
aux.qx.4 <- fit4.qx$mean

res2 <- (dx - (aux.qx.2 * ex))^2 / (aux.qx.2 * ex)
sum(res2, na.rm = T)   ### 38955.08
res3 <- (dx - (aux.qx.3 * ex))^2 / (aux.qx.3 * ex)
sum(res3, na.rm = T)   ### 43483
res4 <- (dx - (aux.qx.4 * ex))^2 / (aux.qx.4 * ex)
sum(res4, na.rm = T)   ### 42494.7

library(BayesMortalityPlus)
fit.mar.rj <- readRDS("rj_fit_mar.RDS")
fit.qx.rj <- fitted(fit.mar.rj)

###nas idades
plot(1:13, 1 - exp(-exp(y_treat[222:234,20])), log = "y", ylim = c(0.001,1))
points(1:13, 1 - exp(-exp(y_treat[1:13, 20])))
points(1:13, 1 - exp(-exp(y_treat[14:26, 20])))
lines(fit4.qx$mean[222:234,20], col = "blue")
lines(fit4.qx$upper[222:234,20], col = "blue", lty = 2)
lines(fit4.qx$lower[222:234,20], col = "blue", lty = 2)

lines(fit3.qx$mean[222:234,20], col = "red")
lines(fit3.qx$upper[222:234,20], col = "red", lty = 2)
lines(fit3.qx$lower[222:234,20], col = "red", lty = 2)

lines(fit.qx.rj$mean[,20], col = "black")
lines(fit.qx.rj$upper[,20], col = "black", lty = 2)
lines(fit.qx.rj$lower[,20], col = "black", lty = 2)

###no ano 
### 20-24
plot(1:42, 1 - exp(-exp(y_treat[222,])), log = "y")
lines(fit4.qx$mean[222,], col = "blue")
lines(fit4.qx$upper[222,], col = "blue", lty = 2)
lines(fit4.qx$lower[222,], col = "blue", lty = 2)

lines(fit3.qx$mean[222,], col = "red")
lines(fit3.qx$upper[222,], col = "red", lty = 2)
lines(fit3.qx$lower[222,], col = "red", lty = 2)

lines(fit.qx.rj$mean[1,], col = "black")
lines(fit.qx.rj$upper[1,], col = "black", lty = 2)
lines(fit.qx.rj$lower[1,], col = "black", lty = 2)
## meio ruim pra todos (melhor pro marginal)

### 35-39
plot(1:42, 1 - exp(-exp(y_treat[225,])), log = "y")
lines(fit4.qx$mean[225,], col = "blue")
lines(fit4.qx$upper[225,], col = "blue", lty = 2)
lines(fit4.qx$lower[225,], col = "blue", lty = 2)

lines(fit3.qx$mean[225,], col = "red")
lines(fit3.qx$upper[225,], col = "red", lty = 2)
lines(fit3.qx$lower[225,], col = "red", lty = 2)

lines(fit.qx.rj$mean[4,], col = "black")
lines(fit.qx.rj$upper[4,], col = "black", lty = 2)
lines(fit.qx.rj$lower[4,], col = "black", lty = 2)

### modelo marginal é obviamente mais flexível que o modelo espacial, nao eh
# comparavel nesse sentido

### sigmas:
1/apply(fit.mar.rj$phiv, 1, median)
median(fit2$sigma_e.chain)
apply(fit3$sigma_e.chain, 2, median) ### espaciais com var_e mais alta mas td bem

1/median(fit.mar.rj$phiw)
median(fit2$sigma_w.chain)
median(fit3$sigma_w.chain) ### menor variancia pra kappa, vai ver melhora forecast

apply(fit2$theta.chain, 2, median)
apply(fit3$theta.chain, 2, median) ## thetas diferentes

apply(fit2$gamma.chain, 2, median)
apply(fit3$gamma.chain, 2, median) ### gammas bastante altos
### testar constraints diferentes em theta e gamma
### implementar predict pra tirar ano de 2021 e tentar prever a frente
### implementar missing com o fit3

colnames(res2) <- 1980:2021
colnames(res3) <- 1980:2021

### RJ
# modelagem univariado, RJ (bom volume de dados) ##33018 -----

dx[which(ind == "X33018.RIO.DE.JANEIRO"),] -> aux.dx ## true
ex[which(ind == "X33018.RIO.DE.JANEIRO"),] -> aux.ex

res_mar <- ((aux.dx - fit.qx.rj$mean*aux.ex)^2)/(fit.qx.rj$mean*aux.ex)
sum(res_mar)   ##12710
sum(res2[which(ind == "X33018.RIO.DE.JANEIRO"),]) ## 18221
sum(res3[which(ind == "X33018.RIO.DE.JANEIRO"),]) ## 23814
sum(res4[which(ind == "X33018.RIO.DE.JANEIRO"),]) ## 22843


## facet_wrap com marginal x conjunto e ver os erros
## dar um jeito de quantificar tb

aux <- rbind(res2[which(ind == "X33018.RIO.DE.JANEIRO"),],
             res3[which(ind == "X33018.RIO.DE.JANEIRO"),], res_mar)

as.data.frame(aux) %>% cbind(age = rep(1:13, 3), id = c(rep(1, 13), rep(2,13), rep(3, 13))) %>%
  pivot_longer(1:42) %>%
  ggplot() +
  geom_tile(aes(x = as.numeric(name), y = age, fill = value)) +
  facet_wrap(~id, ncol = 3, labeller = labeller(id = c("1" = "SBLC",
                                                       "2" = "SBLC sigma_e_age",
                                                       "3" = "BLC"))) +
  scale_fill_gradientn(colors = rev(rainbow(4)))
  #scale_fill_gradient2(low = "white", mid = "yellow", high = "red4", midpoint = 500)




### ST MAR MADALENA
# modelagem univariado, santa maria madalena (pouco volume de dados) ## 33008 ----
library(BayesMortalityPlus)
fit.smm <- readRDS("smm_fit_mar.RDS")
fit.qx.smm <- fitted(fit.smm)

plot(1:13, 1 - exp(-exp(y_treat[92:104,35])), log = "y")
lines(fit4.qx$mean[92:104,35], col = "blue")
lines(fit4.qx$upper[92:104,35], col = "blue", lty = 2)
lines(fit4.qx$lower[92:104,35], col = "blue", lty = 2)

lines(fit3.qx$mean[92:104,35], col = "red")
lines(fit3.qx$upper[92:104,35], col = "red", lty = 2)
lines(fit3.qx$lower[92:104,35], col = "red", lty = 2)

lines(fit.qx.smm$mean[,35], col = "black")
lines(fit.qx.smm$upper[,35], col = "black", lty = 2)
lines(fit.qx.smm$lower[,35], col = "black", lty = 2)

dx[which(ind == "X33008.SANTA.MARIA.MADALENA"),] -> aux.dx ## true
ex[which(ind == "X33008.SANTA.MARIA.MADALENA"),] -> aux.ex

res_mar <- ((aux.dx - fit.qx.smm$mean*aux.ex)^2)/(fit.qx.smm$mean*aux.ex)
sum(res_mar, na.rm = T)  ## 618.8
sum(res2[which(ind == "X33008.SANTA.MARIA.MADALENA"),], na.rm = T) ## 669.8
sum(res3[which(ind == "X33008.SANTA.MARIA.MADALENA"),], na.rm = T) ## 677.9
sum(res4[which(ind == "X33008.SANTA.MARIA.MADALENA"),], na.rm = T) ## 678.7

### modelo com theta = 0 vai melhorar pro espaço fixo, e piorar pro restante


aux <- rbind(res2[which(ind == "X33008.SANTA.MARIA.MADALENA"),],
             res3[which(ind == "X33008.SANTA.MARIA.MADALENA"),], res_mar)


as.data.frame(aux) %>% cbind(age = rep(1:13, 3), id = c(rep(1, 13), rep(2,13), rep(3, 13))) %>%
  pivot_longer(1:42) %>%
  ggplot() +
  geom_tile(aes(x = as.numeric(name), y = age, fill = value)) +
  facet_wrap(~id, ncol = 3, labeller = labeller(id = c("1" = "SBLC",
                                                       "2" = "SBLC sigma_e_age",
                                                       "3" = "BLC"))) +
  scale_fill_gradientn(colors = rev(rainbow(4)))
#### ruim, nao teve ganho de informacao, sigma_e variando sempre pior
#### nova constraint?
