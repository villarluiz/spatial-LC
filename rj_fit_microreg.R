library(sf)
library(spdep)
library(tidyverse)

#rj_mun <- st_read("RJ_Municipios_2022/RJ_Municipios_2022.shp")
rj_micro <- st_read("RJ_Microrregioes_2022/RJ_Microrregioes_2022.shp")
#plot(rj_mun[,c(2,5)])
plot(rj_micro[,c(2,5)])

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

fit2 <- sblc_2(y_treat, ages, t, n, m0, C0, M, W, it = 20000, bn = 10000, thin = 4)
#saveRDS(fit2, "rj_micro_fit2.RDS")
saveRDS(fit2, "rj_micro_fit2_full.RDS")
#### melhorou os parametros essa restricao que fiz pensando no blc.R()

#saveRDS(fit, "rj_micro_fit.RDS")
fit <- readRDS("rj_micro_fit.RDS")

### vis. cadeias

plot.ts(fit$kappa.chain[ ,1])
plot.ts(fit$kappa.chain[ ,10])
plot.ts(fit$kappa.chain[ ,20])
plot.ts(fit$kappa.chain[ ,30])

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

acf(fit$theta.chain[ ,1])
acf(fit$theta.chain[ ,10])
acf(fit$theta.chain[ ,13])
graphics.off()


plot.ts(fit$nu.chain)
acf(fit$nu.chain)
plot.ts(fit$sigma_w.chain)
acf(fit$sigma_w.chain)
plot.ts(fit$sigma_e.chain)
acf(fit$sigma_e.chain)
plot.ts(fit$sigma_t.chain)
acf(fit$sigma_t.chain)

plot.ts(apply(fit$theta.chain, 2, median))
### ok sigma_t ser 5?


plot(apply(fit$kappa.chain, 2, median), type = "l")
lines(apply(fit$kappa.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$kappa.chain, 2, quantile, 0.975), col = 2, lty = 2)

plot(apply(fit$alpha.chain, 2, median), type = "l")
lines(apply(fit$alpha.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$alpha.chain, 2, quantile, 0.975), col = 2, lty = 2)

plot(apply(fit$beta.chain, 2, median), type = "l", ylim = c(0, 0.008))
lines(apply(fit$beta.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$beta.chain, 2, quantile, 0.975), col = 2, lty = 2)
abline(h = 0)

plot(apply(fit$gamma.chain, 2, median), type = "l")
lines(apply(fit$gamma.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$gamma.chain, 2, quantile, 0.975), col = 2, lty = 2)

plot(apply(fit$theta.chain, 2, median), type = "l")
lines(apply(fit$theta.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$theta.chain, 2, quantile, 0.975), col = 2, lty = 2)


#### gerar y_ts
source("fitted_sblc.R")
fit <- readRDS("rj_micro_fit.RDS")

fitted <- fitted_sblc(fit)


## comparar qxs via metrica residual pearson quadrada = (o - e)^2/e, fazer grafico heatmap por reg
library(ggplot2)
aux.qx <- 1 - exp(-exp(fitted$mean))
y_qx <- 1 - exp(-exp(y_treat))
res <- ( (y_qx - aux.qx)^2)/aux.qx
colnames(res) <- 1980:2021

### RJ
# modelagem univariado, RJ (bom volume de dados) ##33018 -----
y <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt"))
  y[,i] <- c(as.numeric(mx$X33018.RIO.DE.JANEIRO))
}
colnames(y) <- as.character(1980:2021)
row.names(y) <- mx[,1]

library(BayesMortalityPlus)
set.seed(14)
fit <- blc(y, M = 10000)
plot(fit)
fitted_mar <- fitted(fit)

saveRDS(fit, "rj_fit_mar.RDS")
# aux.lnmx <- log( -log(1 - fitted_mar$mean) )
# res_mar <- ((y - aux.lnmx)^2)/aux.lnmx

y_treat[which(ind == "X33018.RIO.DE.JANEIRO"),] == y ## true
y2_qx <- 1 - exp(-exp(y))

res_mar <- ((y2_qx - fitted_mar$mean)^2)/fitted_mar$mean


## facet_wrap com marginal x conjunto e ver os erros
## dar um jeito de quantificar tb

rj_fit <- res[which(ind == "X33018.RIO.DE.JANEIRO"),]

aux <- rbind(rj_fit, res_mar)

as.data.frame(aux) %>% cbind(age = rep(1:13, 2), id = c(rep(1, 13), rep(2,13))) %>%
  pivot_longer(1:42) %>%
  ggplot() +
  geom_tile(aes(x = name, y = age, fill = log(value, base = 10))) +
  facet_wrap(~id, ncol = 2, labeller = labeller(id = c("1" = "SBLC",
                                                       "2" = "BLC"))) +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red4", midpoint = -7)




### ST MAR MADALENA
# modelagem univariado, santa maria madalena (pouco volume de dados) ## 33008 ----
y <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt"))
  y[,i] <- c(as.numeric(mx$X33008.SANTA.MARIA.MADALENA))
}
colnames(y) <- as.character(1980:2021)
row.names(y) <- mx[,1]

#### inputacao de dados 
y[1,6] <- mean(y[1,c(5,7)])
y[1,12] <- mean(y[1,c(11,13)])
y[2,12] <- mean(y[2,c(11,14)]); y[2,13] <- mean(y[2,c(12,14)])
y[3,23] <- mean(y[3,c(22,24)])
y[3,23] <- mean(y[3,c(22,24)])
y[1,39] <- mean(y[1,c(38,40)])
y[1,42] <- y[1,41]


set.seed(14)
fit.smm <- blc(y, M = 10000)
plot(fit.smm)

saveRDS(fit.smm, "smm_fit_mar.RDS")
fit.smm <- readRDS("smm_fit_mar.RDS")

sta_fit <- res[which(ind == "X33008.SANTA.MARIA.MADALENA"),]

fitted_mar_smm <- fitted(fit.smm)

y_treat[which(ind == "X33008.SANTA.MARIA.MADALENA"),] == y ## true
y2_qx <- 1 - exp(-exp(y))

res_mar <- ((y2_qx - fitted_mar_smm$mean)^2)/fitted_mar_smm$mean

aux <- rbind(sta_fit, res_mar)

as.data.frame(aux) %>% cbind(age = rep(1:13, 2), id = c(rep(1, 13), rep(2,13))) %>%
  pivot_longer(1:42) %>%
  ggplot() +
  geom_tile(aes(x = name, y = age, fill = log(value, base = 10))) +
  facet_wrap(~id, ncol = 2, labeller = labeller(id = c("1" = "SBLC",
                                                       "2" = "BLC"))) +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red4", midpoint = -7)

#### ruim, nao teve ganho de informacao


ex <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- as.data.frame(t(read.table(paste("Dados_microrregioes/det/ex_det/ex_", as.character(aux), ".txt", sep = ""),sep = ";")))[-c(1,14),20]
  ex[,i] <- as.numeric(mx)
}
colnames(y) <- as.character(1980:2021)


dx <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- as.data.frame(t(read.table(paste("Dados_microrregioes/det/obt_det/obt_", as.character(aux), ".txt", sep = ""),sep = ";")))[-c(1,14),20]
  dx[,i] <- as.numeric(mx)
}
colnames(y) <- as.character(1980:2021)
y <- log(dx/ex)

library(BayesMortalityPlus)
fit.mar.agg <- blc(y)
plot(fit.mar.agg)

median(fit.mar.agg$theta)
median(fit$nu.chain)
