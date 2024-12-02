library(sf)
library(spdep)
library(tidyverse)

#rj_mun <- st_read("RJ_Municipios_2022/RJ_Municipios_2022.shp")
jp <- st_read("jpn_adm_2019_shp/jpn_admbnda_adm1_2019.shp")
#plot(rj_mun[,c(2,5)])
#plot(rj_micro[,c(2,5)])

shape_aisp_nb_val <- st_make_valid(jp[,c(4,7)]) %>%
  mutate(ADM1_EN = str_trim(ADM1_EN)) %>%
  filter(! ADM1_EN %in% c("Hokkaido", "Okinawa"))
nb <- poly2nb(shape_aisp_nb_val, queen=TRUE)

# matriz de vizinhaça
W <- nb2mat(nb, style="B", zero.policy = TRUE)
M <- diag(rowSums(W))

aux.names <- unique(shape_aisp_nb_val$ADM1_EN)
aux.names[12] <- "Hyogo"
aux.names[19] <- "Kochi"
aux.names[29] <- "Oita"
aux.names[10] <- "Gumma"

##subpopulacao de 18 (mesmo do rj)
W.2 <- W[1:18,1:18]
###7,12,17,18 sem vizinhos -> 19, 20, 21, 22
aux.2 <- 1:18
aux.2[c(7,12,17,18)] <- c(19,23,21,22) 
W.2 <- W[aux.2, aux.2]
M.2 <- diag(rowSums(W.2))

aux.names.2 <- aux.names[aux.2]

##fit idades adultas: 20-33
ages <- length(20:33)  ## idades    20-24, 25-29, ..., 80+
t <- 2004-1975+1  ## tempo
n <- 18  ## regioes

it = 2000 #iterations
bn = 1000
thin = 1

m0 <- 0; C0 <- 100

### y_obs
y <- matrix(NA, nrow = ages*n, ncol = t)
for(i in 1:n){
  aux <- aux.names[i]
  mx <- as.data.frame(read.table(paste0("japan_data/", str_trim(as.character(aux)), "_1x1.txt"), skip = 1, header = T))[,-c(3,4)] %>%
    filter(Year %in% 1975:2004, Age %in% 25:40) %>%
    mutate(Total = log(as.numeric(Total))) %>%
    pivot_wider(names_from = Year, values_from = Total)
  y[((1:16) + 16*(i-1)),] <- as.matrix(mx[,-1])
}
ind <- rep(aux.names, each = 100)



### missing data input (interpolando no ano, nao no espaço)
y[!is.finite(y)] <- NA
ind_missing = which(is.na(y), arr.ind = T)
m.inp <- rowMeans(y, na.rm = T)[ind_missing[,1]]

y_treat <- y
y_treat[ind_missing] <- m.inp


### moran
result <- matrix(NA, nrow = ages, ncol = t)
data_w <- nb2listw(nb, style = "W")
for(k in 1:t){
  for(i in 1:ages){
    aux <- y_treat[seq(i, 4500, by = ages), k]
    result[i,k] <- moran.test(aux, data_w)$p.value
  }
}

options(scipen=999)
colnames(result) <- 1975:2004
result.2 <- cbind(as.data.frame(result), x = 0:99) %>%
  pivot_longer(1:t)
plotly::ggplotly(ggplot(result.2) +
                   scale_x_discrete(breaks = scales::pretty_breaks()) +
                   geom_tile(aes(x = name, y = x, fill = value)) +
                   scale_fill_gradient2(low = "red4", mid = "white", high = "blue", midpoint = 0.05) +
                   labs(x = "year", y = "age interval", fill = "p-valor") )

sum(result < 0.05)/(100*30)  ## 34.7% tem correlacao
#### parecido com cenario RJ no indice de Moran



source("sblc_fun.R")
source("sffbs_fun_sigmae.R")
fit3 <- sblc_3(y_treat, ages, t, n, m0, C0, M.2, W.2, it = 2000, bn = 1000, thin = 1)

saveRDS(fit3, "jp_micro_fit3.RDS")
fit <- readRDS("jp_micro_fit3.RDS")
### vis. cadeias

plot.ts(fit$kappa.chain[ ,1])
plot.ts(fit$kappa.chain[ ,10])
plot.ts(fit$kappa.chain[ ,20])
plot.ts(fit$kappa.chain[ ,30])
plot.ts(fit$kappa.chain[ ,42])

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

plot(apply(fit$alpha.chain, 2, median), type = "l")
lines(apply(fit$alpha.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$alpha.chain, 2, quantile, 0.975), col = 2, lty = 2)

plot(apply(fit$beta.chain, 2, median), type = "l")
lines(apply(fit$beta.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$beta.chain, 2, quantile, 0.975), col = 2, lty = 2)
abline(h = 0)

plot(apply(fit$gamma.chain, 2, median), type = "l")
lines(apply(fit$gamma.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$gamma.chain, 2, quantile, 0.975), col = 2, lty = 2)

plot(apply(fit$theta.chain, 2, median), type = "l")
lines(apply(fit$theta.chain, 2, quantile, 0.025), col = 2, lty = 2)
lines(apply(fit$theta.chain, 2, quantile, 0.975), col = 2, lty = 2)

### maior intervalo (interessante)
### deslocado quando comparado com theta_n = 0



#### gerar y_ts
source("fitted_sblc.R")   ###funciona tambem pra sigmae variando na idade
fit.qx <- fitted_sblc(fit)

## comparar qxs via metrica residual pearson quadrada = (o - e)^2/e, fazer grafico heatmap por reg
library(ggplot2)
### y_obs

aux.qx.2 <- fit.qx$mean

res2 <- ( (1 - exp(-exp(y_treat))) - aux.qx.2)^2 / (aux.qx.2)
sum(res2, na.rm = T)   ### 38955.08


library(BayesMortalityPlus)
fit.mar.aichi <- blc(y_treat[1:16,])
fit.qx.aichi <- fitted(fit.mar.aichi)

###nas idades
plot(1:16, 1 - exp(-exp(y_treat[1:16,20])), log = "y")
lines(fit.qx$mean[1:16,20], col = "blue")
lines(fit.qx$upper[1:16,20], col = "blue", lty = 2)
lines(fit.qx$lower[1:16,20], col = "blue", lty = 2)

lines(fit.qx.aichi$mean[,20], col = "red")
lines(fit.qx.aichi$upper[,20], col = "red", lty = 2)
lines(fit.qx.aichi$lower[,20], col = "red", lty = 2)
### ok, marginak ainda melhor


###no ano 
### 25
plot(1:30, 1 - exp(-exp(y_treat[1,])), log = "y")
lines(fit.qx$mean[1,], col = "blue")
lines(fit.qx$upper[1,], col = "blue", lty = 2)
lines(fit.qx$lower[1,], col = "blue", lty = 2)

lines(fit.qx.aichi$mean[1,], col = "red")
lines(fit.qx.aichi$upper[1,], col = "red", lty = 2)
lines(fit.qx.aichi$lower[1,], col = "red", lty = 2)
## melhor pro sblc!!!

### 30
plot(1:30, 1 - exp(-exp(y_treat[6,])), log = "y")
lines(fit.qx$mean[6,], col = "blue")
lines(fit.qx$upper[6,], col = "blue", lty = 2)
lines(fit.qx$lower[6,], col = "blue", lty = 2)

lines(fit.qx.aichi$mean[6,], col = "red")
lines(fit.qx.aichi$upper[6,], col = "red", lty = 2)
lines(fit.qx.aichi$lower[6,], col = "red", lty = 2)

### 35
plot(1:30, 1 - exp(-exp(y_treat[11,])), log = "y")
lines(fit.qx$mean[11,], col = "blue")
lines(fit.qx$upper[11,], col = "blue", lty = 2)
lines(fit.qx$lower[11,], col = "blue", lty = 2)

lines(fit.qx.aichi$mean[11,], col = "red")
lines(fit.qx.aichi$upper[11,], col = "red", lty = 2)
lines(fit.qx.aichi$lower[11,], col = "red", lty = 2)
## melhor pro sblc!!!

### 40
plot(1:30, 1 - exp(-exp(y_treat[16,])), log = "y")
lines(fit.qx$mean[16,], col = "blue")
lines(fit.qx$upper[16,], col = "blue", lty = 2)
lines(fit.qx$lower[16,], col = "blue", lty = 2)

lines(fit.qx.aichi$mean[16,], col = "red")
lines(fit.qx.aichi$upper[16,], col = "red", lty = 2)
lines(fit.qx.aichi$lower[16,], col = "red", lty = 2)
### modelo espacial ficou melhor pros dados nao intervalados, que estranho

### FAZER GRADIENTE


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
