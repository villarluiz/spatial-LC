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
y_treat[92, 42] <- y_treat[92, 41]

source("sffbs_fun_sigmae.R")
fit3 <- readRDS("rj_micro_fit3_full.RDS")

source("fitted_sblc.R")
fit3.qx <- fitted_sblc(fit3)

source("predict_sblc.R")
sblc.pred <- predict.SBLC(fit3, h = 16)
saveRDS(sblc.pred, "rj_micro_fit3_full_pred.RDS")
sblc.pred <- readRDS("rj_micro_fit3_full_pred.RDS")

qts.pred <- apply(sblc.pred$y, 2:3, quantile, 0.5, na.rm = T)
qi.pred <- apply(sblc.pred$y, 2:3, quantile, 0.025, na.rm = T)
qs.pred <- apply(sblc.pred$y, 2:3, quantile, 0.975, na.rm = T)

library(BayesMortalityPlus)
fit.mar.rj <- readRDS("rj_fit_mar.RDS")
rj.pred <- predict(fit.mar.rj, h = 16)
fit.qx.rj <- fitted(fit.mar.rj)

qts.pred.rj <- apply(rj.pred$y, 2:3, quantile, 0.5)
qi.pred.rj <- apply(rj.pred$y, 2:3, quantile, 0.025)
qs.pred.rj <- apply(rj.pred$y, 2:3, quantile, 0.975)

###no ano 
### 20-24
# plot(1:52, c(1 - exp(-exp(y_treat[222,])), rep(NA,10)), log = "y")
# 
# lines(fit3.qx$mean[222,], col = "red")
# lines(fit3.qx$upper[222,], col = "red", lty = 2)
# lines(fit3.qx$lower[222,], col = "red", lty = 2)
# 
# lines(fit.qx.rj$mean[1,], col = "black")
# lines(fit.qx.rj$upper[1,], col = "black", lty = 2)
# lines(fit.qx.rj$lower[1,], col = "black", lty = 2)
# ## meio ruim pra todos (melhor pro marginal)

### 35-39
plot(1:58, c(1 - exp(-exp(y_treat[225,])), rep(NA,16)), log = "y")

lines(1:42, fit3.qx$mean[225,], col = "red")
lines(1:42, fit3.qx$upper[225,], col = "red", lty = 2)
lines(1:42, fit3.qx$lower[225,], col = "red", lty = 2)

lines(43:58, 1 - exp(-exp(qts.pred[ ,225])), col = "red")
lines(43:58, 1 - exp(-exp(qi.pred[ ,225])), col = "red", lty = 2)
lines(43:58, 1 - exp(-exp(qs.pred[ ,225])), col = "red", lty = 2)

lines(fit.qx.rj$mean[4,], col = "black")
lines(fit.qx.rj$upper[4,], col = "black", lty = 2)
lines(fit.qx.rj$lower[4,], col = "black", lty = 2)

lines(43:58, 1 - exp(-exp(qts.pred.rj[ , 4])), col = "black")
lines(43:58, 1 - exp(-exp(qi.pred.rj[ , 4])), col = "black", lty = 2)
lines(43:58, 1 - exp(-exp(qs.pred.rj[ , 4])), col = "black", lty = 2)

### intervalos muito largos, extrapolacao ta a msm coisa, nao ganhou nada


### 45-49
plot(1:58, c(1 - exp(-exp(y_treat[227,])), rep(NA,16)), log = "y")

lines(1:42, fit3.qx$mean[227,], col = "red")
lines(1:42, fit3.qx$upper[227,], col = "red", lty = 2)
lines(1:42, fit3.qx$lower[227,], col = "red", lty = 2)

lines(43:58, 1 - exp(-exp(qts.pred[ ,227])), col = "red")
lines(43:58, 1 - exp(-exp(qi.pred[ ,227])), col = "red", lty = 2)
lines(43:58, 1 - exp(-exp(qs.pred[ ,227])), col = "red", lty = 2)

lines(fit.qx.rj$mean[6,], col = "black")
lines(fit.qx.rj$upper[6,], col = "black", lty = 2)
lines(fit.qx.rj$lower[6,], col = "black", lty = 2)

lines(43:58, 1 - exp(-exp(qts.pred.rj[ , 6])), col = "black")
lines(43:58, 1 - exp(-exp(qi.pred.rj[ , 6])), col = "black", lty = 2)
lines(43:58, 1 - exp(-exp(qs.pred.rj[ , 6])), col = "black", lty = 2)

### intervalos muito largos, extrapolacao ta a msm coisa, nao ganhou nada

