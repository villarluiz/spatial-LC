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
t <- 40  ## tempo (até 2019)
n <- 18   ## regioes

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
y.control <- y

###### SMM missing ----
## holdout de 2015-2019 pra validação
source("sblc_fun.R")
source("sffbs_fun_sigmae.R")
fit <- sblc_missing_2(y[,1:35], ages, (t-5), n, m0, C0, M, W, it = 10000, bn = 5000, thin = 2)
saveRDS(fit, "rj_micro_fit_model1_full.RDS")
fit <- readRDS("rj_micro_fit_model1_full.RDS")



###### ~5% missing   10206 total data ----
10206*0.05  ##510 missing
set.seed(15)
aux <- cbind(sample(1:234, 510, replace=T),
             sample(1:42, 510, replace = T))
y[aux] <- NA

source("sblc_fun.R")
source("sffbs_fun.R")
fit <- sblc_missing(y, ages, t, n, m0, C0, M, W, it = 3000, bn = 1500, thin = 1)
saveRDS(fit, "rj_micro_fit_missing_5percent.RDS")
fit <- readRDS("rj_micro_fit_missing_5percent.RDS")

k = 3
x11()
par(mfrow=c(3,6))

## input da media temporal na idade 
ind_missing = which(is.na(y), arr.ind = T)
m.inp <- rowMeans(y, na.rm = T)[ind_missing[,1]]

for(i in (1:18) + 18*k){
  plot.ts(fit$input.chain[,i])
  abline(h = y.control[ind_missing][i], col = "red")
  abline(h = m.inp[i], col = "blue")
}

## media dos inputs
mu.y <- apply(fit$input.chain, 2, median)

for(i in (1:18) + 18*k){
  plot(0:1, rep(mu.y[i],2), type = "l")
  abline(h = y.control[ind_missing][i], col = "red")
  abline(h = aux[i], col = "blue")
}

graphics.off()

zoo1 <- zoo::na.approx(y) ###interpola os proximos
zoo1.2 <- zoo::na.spline(y) ###interpola via spline
zoo2 <- zoo::na.aggregate(y) ###substitui pela media da linha (parece com rowMeans)
zoo3 <- t(zoo::na.locf(t(y))) ### repete a ultima obs (naive)

res.sblc <- (exp(y.control[ind_missing]) - exp(mu.y))^2/exp(mu.y)
res.age.inp <- (exp(y.control[ind_missing]) - exp(m.inp))^2/exp(m.inp)

res.zoo1 <- (exp(y.control[ind_missing]) - exp(zoo1[ind_missing]))^2/exp(zoo1[ind_missing])
res.zoo1.2 <- (exp(y.control[ind_missing]) - exp(zoo1.2[ind_missing]))^2/exp(zoo1.2[ind_missing])
res.zoo2 <- (exp(y.control[ind_missing]) - exp(zoo2[ind_missing]))^2/exp(zoo2[ind_missing])
res.zoo3 <- (exp(y.control[ind_missing]) - exp(zoo3[ind_missing]))^2/exp(zoo3[ind_missing])

sum(res.sblc < res.age.inp, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das inputacoes do rowMeans

sum(res.sblc < res.zoo1, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das interpolacoes dos proxs

sum(res.sblc < res.zoo1.2, na.rm = T)/nrow(ind_missing)
### ~27% dos casos ele tem menor residuo de pearson das splines

sum(res.sblc < res.zoo2, na.rm = T)/nrow(ind_missing)
### ~41% dos casos ele tem menor residuo de pearson das inputacoes do zoo (media)

sum(res.sblc < res.zoo3, na.rm = T)/nrow(ind_missing)
### ~15% dos casos ele tem menor residuo de pearson dos naive no tempo



###### 30% missing em um região ----
546*0.3  ##~162 missing
set.seed(15)
aux <- cbind(sample(1:13, 162, replace=T),
             sample(1:42, 162, replace = T))
y[aux] <- NA

source("sblc_fun.R")
source("sffbs_fun_sigmae.R")
fit <- sblc_missing_2(y, ages, t, n, m0, C0, M, W, it = 2000, bn = 1000, thin = 1)
saveRDS(fit, "rj_micro_fit_missing_30ita.RDS")
fit <- readRDS("rj_micro_fit_missing_30ita.RDS")

source("blc_missing_mar.R")
source("kd_filter.R")
source("kd_smoother.R")
library(MASS)
fit.mar <- blc_missing(y[1:13,], M = 2000, bn = 1000, thin = 1)

k = 3
x11()
par(mfrow=c(3,6))

## input da media temporal na idade 
ind_missing = which(is.na(y), arr.ind = T)
m.inp <- rowMeans(y, na.rm = T)[ind_missing[,1]]

for(i in (1:18) + 18*k){
  plot.ts(fit$input.chain[,i])
  abline(h = y.control[ind_missing][i], col = "red")
  abline(h = m.inp[i], col = "blue")
}

## media dos inputs
mu.y <- apply(fit$input.chain, 2, median)

zoo1 <- zoo::na.approx(y) ###interpola os proximos
zoo1.2 <- zoo::na.spline(y) ###interpola via spline
zoo2 <- zoo::na.aggregate(y) ###substitui pela media da linha (parece com rowMeans)
zoo3 <- t(zoo::na.locf(t(y))) ### repete a ultima obs (naive)

res.sblc <- (exp(y.control[ind_missing]) - exp(mu.y))^2/exp(mu.y)
res.age.inp <- (exp(y.control[ind_missing]) - exp(m.inp))^2/exp(m.inp)

res.zoo1 <- (exp(y.control[ind_missing]) - exp(zoo1[ind_missing]))^2/exp(zoo1[ind_missing])
res.zoo1.2 <- (exp(y.control[ind_missing]) - exp(zoo1.2[ind_missing]))^2/exp(zoo1.2[ind_missing])
res.zoo2 <- (exp(y.control[ind_missing]) - exp(zoo2[ind_missing]))^2/exp(zoo2[ind_missing])
res.zoo3 <- (exp(y.control[ind_missing]) - exp(zoo3[ind_missing]))^2/exp(zoo3[ind_missing])

sum(res.sblc < res.age.inp, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das inputacoes do rowMeans

sum(res.sblc < res.zoo1, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das interpolacoes dos proxs

sum(res.sblc < res.zoo1.2, na.rm = T)/nrow(ind_missing)
### ~27% dos casos ele tem menor residuo de pearson das splines

sum(res.sblc < res.zoo2, na.rm = T)/nrow(ind_missing)
### ~41% dos casos ele tem menor residuo de pearson das inputacoes do zoo (media)

sum(res.sblc < res.zoo3, na.rm = T)/nrow(ind_missing)
### ~15% dos casos ele tem menor residuo de pearson dos naive no tempo

source("fitted_sblc.R")
qx.sblc <- fitted_sblc(fit)
library(BayesMortalityPlus)
qx.mar <- fitted(fit.mar)
plot(1 - exp(-exp(y.control[1, ])), log = "y")
points(1 - exp(-exp(y[1,])), pch = 16)
lines(qx.sblc$mean[1,], col  ='red')
lines(qx.sblc$upper[1,], col  ='red', lty = 2)
lines(qx.sblc$lower[1,], col  ='red', lty = 2)

lines(qx.mar$mean[1,], col  ='blue')
lines(qx.mar$upper[1,], col  ='blue', lty = 2)
lines(qx.mar$lower[1,], col  ='blue', lty = 2)



plot(1 - exp(-exp(y.control[3, ])), log = "y")
points(1 - exp(-exp(y[3,])), pch = 16)
lines(qx.sblc$mean[3,], col  ='red')
lines(qx.sblc$upper[3,], col  ='red', lty = 2)
lines(qx.sblc$lower[3,], col  ='red', lty = 2)

lines(qx.mar$mean[3,], col  ='blue')
lines(qx.mar$upper[3,], col  ='blue', lty = 2)
lines(qx.mar$lower[3,], col  ='blue', lty = 2)
#### avaliar

###### 50% missing em um região ----
546*0.5  ##~162 missing
set.seed(15)
aux <- cbind(sample(1:13, 273, replace=T),
             sample(1:42, 273, replace = T))
y[aux] <- NA

source("sblc_fun.R")
source("sffbs_fun_sigmae.R")
fit <- sblc_missing_2(y, ages, t, n, m0, C0, M, W, it = 3000, bn = 1500, thin = 1)
saveRDS(fit, "rj_micro_fit_missing_50ita.RDS")
fit <- readRDS("rj_micro_fit_missing_50ita.RDS")

source("blc_missing_mar.R")
source("kd_filter.R")
source("kd_smoother.R")
library(MASS)
fit.mar <- blc_missing(y[1:13,], M = 2000, bn = 1000, thin = 1)

k = 3
x11()
par(mfrow=c(3,6))

## input da media temporal na idade 
ind_missing = which(is.na(y), arr.ind = T)
m.inp <- rowMeans(y, na.rm = T)[ind_missing[,1]]

for(i in (1:18) + 18*k){
  plot.ts(fit$input.chain[,i])
  abline(h = y.control[ind_missing][i], col = "red")
  abline(h = m.inp[i], col = "blue")
}

## media dos inputs
mu.y <- apply(fit$input.chain, 2, median)

zoo1 <- zoo::na.approx(y) ###interpola os proximos
zoo1.2 <- zoo::na.spline(y) ###interpola via spline
zoo2 <- zoo::na.aggregate(y) ###substitui pela media da linha (parece com rowMeans)
zoo3 <- t(zoo::na.locf(t(y))) ### repete a ultima obs (naive)

res.sblc <- (exp(y.control[ind_missing]) - exp(mu.y))^2/exp(mu.y)
res.age.inp <- (exp(y.control[ind_missing]) - exp(m.inp))^2/exp(m.inp)

res.zoo1 <- (exp(y.control[ind_missing]) - exp(zoo1[ind_missing]))^2/exp(zoo1[ind_missing])
res.zoo1.2 <- (exp(y.control[ind_missing]) - exp(zoo1.2[ind_missing]))^2/exp(zoo1.2[ind_missing])
res.zoo2 <- (exp(y.control[ind_missing]) - exp(zoo2[ind_missing]))^2/exp(zoo2[ind_missing])
res.zoo3 <- (exp(y.control[ind_missing]) - exp(zoo3[ind_missing]))^2/exp(zoo3[ind_missing])

sum(res.sblc < res.age.inp, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das inputacoes do rowMeans

sum(res.sblc < res.zoo1, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das interpolacoes dos proxs

sum(res.sblc < res.zoo1.2, na.rm = T)/nrow(ind_missing)
### ~27% dos casos ele tem menor residuo de pearson das splines

sum(res.sblc < res.zoo2, na.rm = T)/nrow(ind_missing)
### ~41% dos casos ele tem menor residuo de pearson das inputacoes do zoo (media)

sum(res.sblc < res.zoo3, na.rm = T)/nrow(ind_missing)
### ~15% dos casos ele tem menor residuo de pearson dos naive no tempo

source("fitted_sblc.R")
qx.sblc <- fitted_sblc(fit)
library(BayesMortalityPlus)
qx.mar <- fitted(fit.mar)
plot(1 - exp(-exp(y.control[1, ])), log = "y")
points(1 - exp(-exp(y[1,])), pch = 16)
lines(qx.sblc$mean[1,], col  ='red')
lines(qx.sblc$upper[1,], col  ='red', lty = 2)
lines(qx.sblc$lower[1,], col  ='red', lty = 2)

lines(qx.mar$mean[1,], col  ='blue')
lines(qx.mar$upper[1,], col  ='blue', lty = 2)
lines(qx.mar$lower[1,], col  ='blue', lty = 2)



plot(1 - exp(-exp(y.control[3, ])), log = "y")
points(1 - exp(-exp(y[3,])), pch = 16)
lines(qx.sblc$mean[3,], col  ='red')
lines(qx.sblc$upper[3,], col  ='red', lty = 2)
lines(qx.sblc$lower[3,], col  ='red', lty = 2)

lines(qx.mar$mean[3,], col  ='blue')
lines(qx.mar$upper[3,], col  ='blue', lty = 2)
lines(qx.mar$lower[3,], col  ='blue', lty = 2)

plot(1 - exp(-exp(y.control[8, ])), log = "y")
points(1 - exp(-exp(y[8,])), pch = 16)
lines(qx.sblc$mean[8,], col  ='red')
lines(qx.sblc$upper[8,], col  ='red', lty = 2)
lines(qx.sblc$lower[8,], col  ='red', lty = 2)

lines(qx.mar$mean[8,], col  ='blue')
lines(qx.mar$upper[8,], col  ='blue', lty = 2)
lines(qx.mar$lower[8,], col  ='blue', lty = 2)


plot(1 - exp(-exp(y.control[12, ])), log = "y")
points(1 - exp(-exp(y[12,])), pch = 16)
lines(qx.sblc$mean[12,], col  ='red')
lines(qx.sblc$upper[12,], col  ='red', lty = 2)
lines(qx.sblc$lower[12,], col  ='red', lty = 2)

lines(qx.mar$mean[12,], col  ='blue')
lines(qx.mar$upper[12,], col  ='blue', lty = 2)
lines(qx.mar$lower[12,], col  ='blue', lty = 2)
## RUIM


###### 2020, 21 e 22 missing ----
y = y.control
y[1:13, 40:42] <- NA 

source("sblc_fun.R")
source("sffbs_fun_sigmae.R")
fit <- sblc_missing_2(y, ages, t, n, m0, C0, M, W, it = 3000, bn = 1500, thin = 1)
saveRDS(fit, "rj_micro_fit_missing_3yearsita.RDS")
fit <- readRDS("rj_micro_fit_missing_3yearsita.RDS")

source("blc_missing_mar.R")
source("kd_filter.R")
source("kd_smoother.R")
library(MASS)
fit.mar <- blc_missing(y[1:13,], M = 2000, bn = 1000, thin = 1)

k = 3
x11()
par(mfrow=c(3,6))

## input da media temporal na idade 
ind_missing = which(is.na(y), arr.ind = T)
m.inp <- rowMeans(y, na.rm = T)[ind_missing[,1]]

for(i in (1:18) + 18*k){
  plot.ts(fit$input.chain[,i])
  abline(h = y.control[ind_missing][i], col = "red")
  abline(h = m.inp[i], col = "blue")
}

## media dos inputs
mu.y <- apply(fit$input.chain, 2, median)

zoo1 <- zoo::na.approx(y) ###interpola os proximos
zoo1.2 <- zoo::na.spline(y) ###interpola via spline
zoo2 <- zoo::na.aggregate(y) ###substitui pela media da linha (parece com rowMeans)
zoo3 <- t(zoo::na.locf(t(y))) ### repete a ultima obs (naive)

res.sblc <- (exp(y.control[ind_missing]) - exp(mu.y))^2/exp(mu.y)
res.age.inp <- (exp(y.control[ind_missing]) - exp(m.inp))^2/exp(m.inp)

res.zoo1 <- (exp(y.control[ind_missing]) - exp(zoo1[ind_missing]))^2/exp(zoo1[ind_missing])
res.zoo1.2 <- (exp(y.control[ind_missing]) - exp(zoo1.2[ind_missing]))^2/exp(zoo1.2[ind_missing])
res.zoo2 <- (exp(y.control[ind_missing]) - exp(zoo2[ind_missing]))^2/exp(zoo2[ind_missing])
res.zoo3 <- (exp(y.control[ind_missing]) - exp(zoo3[ind_missing]))^2/exp(zoo3[ind_missing])

sum(res.sblc < res.age.inp, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das inputacoes do rowMeans

sum(res.sblc < res.zoo1, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das interpolacoes dos proxs

sum(res.sblc < res.zoo1.2, na.rm = T)/nrow(ind_missing)
### ~27% dos casos ele tem menor residuo de pearson das splines

sum(res.sblc < res.zoo2, na.rm = T)/nrow(ind_missing)
### ~41% dos casos ele tem menor residuo de pearson das inputacoes do zoo (media)

sum(res.sblc < res.zoo3, na.rm = T)/nrow(ind_missing)
### ~15% dos casos ele tem menor residuo de pearson dos naive no tempo

source("fitted_sblc.R")
qx.sblc <- fitted_sblc(fit)
library(BayesMortalityPlus)
qx.mar <- fitted(fit.mar)
plot(1 - exp(-exp(y.control[1, ])), log = "y")
points(1 - exp(-exp(y[1,])), pch = 16)
lines(qx.sblc$mean[1,], col  ='red')
lines(qx.sblc$upper[1,], col  ='red', lty = 2)
lines(qx.sblc$lower[1,], col  ='red', lty = 2)

lines(qx.mar$mean[1,], col  ='blue')
lines(qx.mar$upper[1,], col  ='blue', lty = 2)
lines(qx.mar$lower[1,], col  ='blue', lty = 2)



plot(1 - exp(-exp(y.control[3, ])), log = "y")
points(1 - exp(-exp(y[3,])), pch = 16)
lines(qx.sblc$mean[3,], col  ='red')
lines(qx.sblc$upper[3,], col  ='red', lty = 2)
lines(qx.sblc$lower[3,], col  ='red', lty = 2)

lines(qx.mar$mean[3,], col  ='blue')
lines(qx.mar$upper[3,], col  ='blue', lty = 2)
lines(qx.mar$lower[3,], col  ='blue', lty = 2)

plot(1 - exp(-exp(y.control[8, ])), log = "y")
points(1 - exp(-exp(y[8,])), pch = 16)
lines(qx.sblc$mean[8,], col  ='red')
lines(qx.sblc$upper[8,], col  ='red', lty = 2)
lines(qx.sblc$lower[8,], col  ='red', lty = 2)

lines(qx.mar$mean[8,], col  ='blue')
lines(qx.mar$upper[8,], col  ='blue', lty = 2)
lines(qx.mar$lower[8,], col  ='blue', lty = 2)


plot(1 - exp(-exp(y.control[12, ])), log = "y")
points(1 - exp(-exp(y[12,])), pch = 16)
lines(qx.sblc$mean[12,], col  ='red')
lines(qx.sblc$upper[12,], col  ='red', lty = 2)
lines(qx.sblc$lower[12,], col  ='red', lty = 2)

lines(qx.mar$mean[12,], col  ='blue')
lines(qx.mar$upper[12,], col  ='blue', lty = 2)
lines(qx.mar$lower[12,], col  ='blue', lty = 2)
## SBLC conseguiu prever a subida dos anos de covid 
## por causa da informacao emprestada 

## RESULTADO BOM ?????



###### 70% missing em um região ----
546*0.7  ##~382 missing
set.seed(15)
aux <- cbind(sample(1:13, 382, replace=T),
             sample(1:42, 382, replace = T))
y[aux] <- NA

source("sblc_fun.R")
source("sffbs_fun_sigmae.R")
fit <- sblc_missing_2(y, ages, t, n, m0, C0, M, W, it = 10000, bn = 5000, thin = 2)
saveRDS(fit, "rj_micro_fit_missing_70ita.RDS")
fit <- readRDS("rj_micro_fit_missing_70ita.RDS")

source("blc_missing_mar.R")
source("kd_filter.R")
source("kd_smoother.R")
library(MASS)
fit.mar <- blc_missing(y[1:13,], M = 2000, bn = 1000, thin = 1)

# k = 3
# x11()
# par(mfrow=c(3,6))

## input da media temporal na idade 
ind_missing = which(is.na(y), arr.ind = T)
m.inp <- rowMeans(y, na.rm = T)[ind_missing[,1]]

# for(i in (1:18) + 18*k){
#   plot.ts(fit$input.chain[,i])
#   abline(h = y.control[ind_missing][i], col = "red")
#   abline(h = m.inp[i], col = "blue")
# }

## media dos inputs
mu.y <- apply(fit$input.chain, 2, median)

zoo1 <- zoo::na.approx(y) ###interpola os proximos
zoo1.2 <- zoo::na.spline(y) ###interpola via spline
zoo2 <- zoo::na.aggregate(y) ###substitui pela media da linha (parece com rowMeans)
zoo3 <- t(zoo::na.locf(t(y))) ### repete a ultima obs (naive)

res.sblc <- (exp(y.control[ind_missing]) - exp(mu.y))^2/exp(mu.y)
res.age.inp <- (exp(y.control[ind_missing]) - exp(m.inp))^2/exp(m.inp)

res.zoo1 <- (exp(y.control[ind_missing]) - exp(zoo1[ind_missing]))^2/exp(zoo1[ind_missing])
res.zoo1.2 <- (exp(y.control[ind_missing]) - exp(zoo1.2[ind_missing]))^2/exp(zoo1.2[ind_missing])
res.zoo2 <- (exp(y.control[ind_missing]) - exp(zoo2[ind_missing]))^2/exp(zoo2[ind_missing])
res.zoo3 <- (exp(y.control[ind_missing]) - exp(zoo3[ind_missing]))^2/exp(zoo3[ind_missing])

sum(res.sblc < res.age.inp, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das inputacoes do rowMeans

sum(res.sblc < res.zoo1, na.rm = T)/nrow(ind_missing)
### ~20% dos casos ele tem menor residuo de pearson das interpolacoes dos proxs

sum(res.sblc < res.zoo1.2, na.rm = T)/nrow(ind_missing)
### ~27% dos casos ele tem menor residuo de pearson das splines

sum(res.sblc < res.zoo2, na.rm = T)/nrow(ind_missing)
### ~41% dos casos ele tem menor residuo de pearson das inputacoes do zoo (media)

sum(res.sblc < res.zoo3, na.rm = T)/nrow(ind_missing)
### ~15% dos casos ele tem menor residuo de pearson dos naive no tempo

source("fitted_sblc.R")
qx.sblc <- fitted_sblc(fit)
library(BayesMortalityPlus)
qx.mar <- fitted(fit.mar)
plot(1 - exp(-exp(y.control[1, ])), log = "y")
points(1 - exp(-exp(y[1,])), pch = 16)
lines(qx.sblc$mean[1,], col  ='red')
lines(qx.sblc$upper[1,], col  ='red', lty = 2)
lines(qx.sblc$lower[1,], col  ='red', lty = 2)

lines(qx.mar$mean[1,], col  ='blue')
lines(qx.mar$upper[1,], col  ='blue', lty = 2)
lines(qx.mar$lower[1,], col  ='blue', lty = 2)



plot(1 - exp(-exp(y.control[3, ])), log = "y")
points(1 - exp(-exp(y[3,])), pch = 16)
lines(qx.sblc$mean[3,], col  ='red')
lines(qx.sblc$upper[3,], col  ='red', lty = 2)
lines(qx.sblc$lower[3,], col  ='red', lty = 2)

lines(qx.mar$mean[3,], col  ='blue')
lines(qx.mar$upper[3,], col  ='blue', lty = 2)
lines(qx.mar$lower[3,], col  ='blue', lty = 2)

plot(1 - exp(-exp(y.control[8, ])), log = "y")
points(1 - exp(-exp(y[8,])), pch = 16)
lines(qx.sblc$mean[8,], col  ='red')
lines(qx.sblc$upper[8,], col  ='red', lty = 2)
lines(qx.sblc$lower[8,], col  ='red', lty = 2)

lines(qx.mar$mean[8,], col  ='blue')
lines(qx.mar$upper[8,], col  ='blue', lty = 2)
lines(qx.mar$lower[8,], col  ='blue', lty = 2)


plot(1 - exp(-exp(y.control[12, ])), log = "y")
points(1 - exp(-exp(y[12,])), pch = 16)
lines(qx.sblc$mean[12,], col  ='red')
lines(qx.sblc$upper[12,], col  ='red', lty = 2)
lines(qx.sblc$lower[12,], col  ='red', lty = 2)

lines(qx.mar$mean[12,], col  ='blue')
lines(qx.mar$upper[12,], col  ='blue', lty = 2)
lines(qx.mar$lower[12,], col  ='blue', lty = 2)
## RUIM

