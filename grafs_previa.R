library(sf)
library(spdep)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(patchwork)

ages <- 13  ## idades    20-24, 25-29, ..., 80+
t <- 40  ## tempo (até 2019)
n <- 18


#### RIO DE JANEIRO ----
#### grafico spatial
rj_micro <- st_read("RJ_Microrregioes_2022/RJ_Microrregioes_2022.shp")
pdf("rj_spatial.pdf", width=7, height=3)
ggplot(rj_micro) + geom_sf(aes(fill = NM_MICRO)) + theme_bw() + labs(fill = "") +
  guides(fill=guide_legend(ncol=2))
graphics.off()

y <- matrix(NA, nrow = ages*n, ncol = t)
for(i in 1:t){
  aux <- i + 1979
  mx <- as.data.frame(read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt")))[,-1] %>%
    stack()
  y[,i] <- mx$values
}
ind <- as.character(mx$ind)

### heatmap 25-29 em 1980-2000-2019
mx <- read.table("Dados_microrregioes/det/logmx_det/logmx_1980.txt")
mx2 <- read.table("Dados_microrregioes/det/logmx_det/logmx_2000.txt")
mx3 <- read.table("Dados_microrregioes/det/logmx_det/logmx_2019.txt")

data.aux <- cbind(rj_micro, t(mx[2,-1]), t(mx2[2,-1]), t(mx3[2,-1]))
data.aux.rang <- cbind(t(mx[2,-1]), t(mx2[2,-1]), t(mx3[2,-1]))

library(ggplot2)
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X2)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("1980") -> plt1
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X2.1)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("2000") -> plt2
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X2.2)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("2019") -> plt3

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

library(patchwork)
plt1 + plt2 + plt3 +
  plot_annotation(title = "Mortality rates for 25-29 age group, RJ", theme = theme(plot.title = element_text(color="black", size=16))) +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        low = "yellow", high = "red") & 
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14)) -> pt
pdf("rj_heatmap_2529.pdf", width=12, height=3)
pt
graphics.off()

## 50-54
data.aux <- cbind(rj_micro, t(mx[7,-1]), t(mx2[7,-1]), t(mx3[7,-1]))
data.aux.rang <- cbind(t(mx[7,-1]), t(mx2[7,-1]), t(mx3[7,-1]))

library(ggplot2)
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X7)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("1980") -> plt4
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X7.1)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("2000") -> plt5
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X7.2)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("2019") -> plt6

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

library(patchwork)
plt4 + plt5 + plt6 +
  plot_annotation(title = "Mortality rates for 50-54 age group, RJ", theme = theme(plot.title = element_text(color="black", size=16))) +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        high = "#00695C", low = "#8BC34A") & 
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14)) -> pt2
pdf("rj_heatmap_5054.pdf", width=12, height=3)
pt2
graphics.off()


## 80+
data.aux <- cbind(rj_micro, t(mx[13,-1]), t(mx2[13,-1]), t(mx3[13,-1]))
data.aux.rang <- cbind(t(mx[13,-1]), t(mx2[13,-1]), t(mx3[13,-1]))

ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X13)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("1980") -> plt7
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X13.1)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("2000") -> plt8
ggplot(data = data.aux) + geom_sf(aes(fill = exp(as.numeric(X13.2)))) + labs(fill = "mx") +
  theme_bw() +
  ggtitle("2019") -> plt9

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

plt7 + plt8 + plt9 +
  plot_annotation(title = "Mortality rates for 80+ age group, RJ", theme = theme(plot.title = element_text(color="black", size=16))) +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        low = "#26C6DA", high = "#9C27B0") & 
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14)) -> pt3
pdf("rj_heatmap_80.pdf", width=12, height=3)
pt3
graphics.off()


### parametros
fit.rj <- readRDS("rj_micro_fit_model1_full.RDS")
# alpha
alp <- t(apply(fit.rj$alpha.chain, 2, quantile, c(0.5, 0.025, 0.975)))
colnames(alp) <- c("fit", "lower", "upper")
ggplot(data = alp) + geom_line(aes(x = 1:13, y = fit)) + theme_bw() +
  geom_ribbon(aes(x = 1:13, ymax = upper, ymin = lower), fill = "grey40", alpha = .4) +
  labs(y = expression(alpha[x]), x = "Age group", fill = "") +
  scale_x_continuous(breaks = c(1, 5, 9, 13), labels = c("20-24", "40-44", "60-64" ,"80+")) -> pt1

# beta
alp <- t(apply(fit.rj$beta.chain, 2, quantile, c(0.5, 0.025, 0.975)))
colnames(alp) <- c("fit", "lower", "upper")
ggplot(data = alp) + geom_line(aes(x = 1:13, y = fit)) + theme_bw() +
  geom_ribbon(aes(x = 1:13, ymax = upper, ymin = lower), fill = "grey40", alpha = .4) +
  labs(y = expression(beta[x]), x = "Age group", fill = "") +
  scale_y_continuous(limits = c(0, 0.022)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 5, 9, 13), labels = c("20-24", "40-44", "60-64" ,"80+")) -> pt2

# gamma
alp <- t(apply(fit.rj$gamma.chain, 2, quantile, c(0.5, 0.025, 0.975)))
colnames(alp) <- c("fit", "lower", "upper")
ggplot(data = alp) + geom_line(aes(x = 1:13, y = fit)) + theme_bw() +
  geom_ribbon(aes(x = 1:13, ymax = upper, ymin = lower), fill = "grey40", alpha = .4) +
  labs(y = expression(gamma[x]), x = "Age group", fill = "") +
  scale_x_continuous(breaks = c(1, 5, 9, 13), labels = c("20-24", "40-44", "60-64" ,"80+")) -> pt3

# kappa
alp <- t(apply(fit.rj$kappa.chain, 2, quantile, c(0.5, 0.025, 0.975)))
colnames(alp) <- c("fit", "lower", "upper")
ggplot(data = alp) + geom_line(aes(x = 1980:2014, y = fit)) + theme_bw() +
  geom_ribbon(aes(x = 1980:2014, ymax = upper, ymin = lower), fill = "grey40", alpha = .4) +
  labs(y = expression(kappa[t]), x = "Year", fill = "") -> pt4

# theta
alp <- apply(fit.rj$theta.chain, 2, median) ##menos diminui mort e positivo aumenta
data.aux <- cbind(rj_micro[,c(2:5)], theta = alp)
ggplot(data.aux) + geom_sf(aes(fill = theta)) + theme_bw() + labs(fill = "") +
  guides(fill = guide_colorbar(title = expression(theta[s]))) +
  scale_fill_gradient2(low = "blue", midpoint = 0, high = "red") -> pt5

pdf("rj_parameters.pdf", width=12, height=4)
design <- "125
           345"
pt1 + pt2 + pt3 + pt4 + free(pt5, side = "tb") + plot_layout(guides = "collect", design = design) &
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14))
graphics.off()

### comp marginal
source("fitted_sblc.R")
qx.sblc.rj <- fitted_sblc(fit.rj)
source("predict_sblc.R")
#qx.sblc.pred <- predict.SBLC(fit.rj, h = 5)
#saveRDS(qx.sblc.pred, "rj_micro_qxpred.RDS")
qx.sblc.pred <- readRDS("rj_micro_qxpred.RDS")
qx.pred <- apply(qx.sblc.pred$y, 2:3, quantile, c(0.5, 0.025, 0.975))


library(BayesMortalityPlus)
fit.mar.rj <- blc(y[222:234, 1:35])
#plot(fit.mar.rj)
qx.m.rj <- fitted(fit.mar.rj)
qx.pred.rj <- fitted(predict(fit.mar.rj, h = 5))

fit.mar.ita <- blc(y[1:13, 1:35])
#plot(fit.mar.ita)
qx.m.ita <- fitted(fit.mar.ita)
qx.pred.ita <- fitted(predict(fit.mar.ita, h = 5))

source("blc_missing_mar.R")
source("kd_filter.R")
source("kd_smoother.R")
library(MASS)
fit.mar.smm <- blc_missing(y[92:104, 1:35])  ### tratamento de missing
plot(fit.mar.smm)
qx.m.smm <- fitted(fit.mar.smm)
source("predict_blc_missing.R")
qx.pred.smm <- fitted(predict_missing_BLC(fit.mar.smm, h = 5))  ###tem que fazer o predict missing
acct <- function(x){log(-log(1 - x))}
# 25-29 RJ
#pdf("rj_rj_qxs.pdf", width=12, height=7)
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.rj$mean[2,]), lmi = acct(qx.m.rj$lower[2,]),
                     lmu = acct(qx.m.rj$upper[2,]), lmx.s = qx.sblc.rj$mean[223,], lmi.s = qx.sblc.rj$lower[223,],
                     lmu.s = qx.sblc.rj$upper[223,], y = y[223,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.rj$mean[2,]), lmi = acct(qx.pred.rj$lower[2,]),
                          lmu = acct(qx.pred.rj$upper[2,]), lmx.s = qx.pred[1, ,223], lmi.s = qx.pred[2, ,223],
                          lmu.s = qx.pred[3, ,223], y = y[223,36:40])

ggplot() + theme_bw() + coord_cartesian(ylim = c(-6.6, -5.7)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "", title = "25-29") -> pt1


# 50-54
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.rj$mean[7,]), lmi = acct(qx.m.rj$lower[7,]),
                     lmu = acct(qx.m.rj$upper[7,]), lmx.s = qx.sblc.rj$mean[228,], lmi.s = qx.sblc.rj$lower[228,],
                     lmu.s = qx.sblc.rj$upper[228,], y = y[228,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.rj$mean[7,]), lmi = acct(qx.pred.rj$lower[7,]),
                          lmu = acct(qx.pred.rj$upper[7,]), lmx.s = qx.pred[1, ,228], lmi.s = qx.pred[2, ,228],
                          lmu.s = qx.pred[3, ,228], y = y[228,36:40])

ggplot() + theme_bw() + coord_cartesian(ylim = c(-5.15, -4.38)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "", title = "50-54") -> pt2

# 80+
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.rj$mean[13,]), lmi = acct(qx.m.rj$lower[13,]),
                     lmu = acct(qx.m.rj$upper[13,]), lmx.s = qx.sblc.rj$mean[234,], lmi.s = qx.sblc.rj$lower[234,],
                     lmu.s = qx.sblc.rj$upper[234,], y = y[234,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.rj$mean[13,]), lmi = acct(qx.pred.rj$lower[13,]),
                          lmu = acct(qx.pred.rj$upper[13,]), lmx.s = qx.pred[1, ,234], lmi.s = qx.pred[2, ,234],
                          lmu.s = qx.pred[3, ,234], y = y[234,36:40])

ggplot() + theme_bw() + coord_cartesian(ylim = c(-2.4, -1.85)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "", title = "80+") -> pt3

# pt1 + pt2 + pt3 +
#   plot_annotation(title = "Rio de Janeiro", theme = theme(plot.title = element_text(color="black", size=16))) +
#   plot_layout(guides = "collect") & 
#   theme(legend.text=element_text(color="black", size=14),
#         legend.title = element_text(color="black", size=14),
#         axis.title = element_text(color = "black", size = 14)) 
# pdf("rj_qx.pdf", width=12, height=3)
# pt
# graphics.off()


### ITA
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.ita$mean[2,]), lmi = acct(qx.m.ita$lower[2,]),
                     lmu = acct(qx.m.ita$upper[2,]), lmx.s = qx.sblc.rj$mean[2,], lmi.s = qx.sblc.rj$lower[2,],
                     lmu.s = qx.sblc.rj$upper[2,], y = y[2,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.ita$mean[2,]), lmi = acct(qx.pred.ita$lower[2,]),
                          lmu = acct(qx.pred.ita$upper[2,]), lmx.s = qx.pred[1, ,2], lmi.s = qx.pred[2, ,2],
                          lmu.s = qx.pred[3, ,2], y = y[2,36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-6.6, -5.7)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt4


# 50-54
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.ita$mean[7,]), lmi = acct(qx.m.ita$lower[7,]),
                     lmu = acct(qx.m.ita$upper[7,]), lmx.s = qx.sblc.rj$mean[7,], lmi.s = qx.sblc.rj$lower[7,],
                     lmu.s = qx.sblc.rj$upper[7,], y = y[7,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.ita$mean[7,]), lmi = acct(qx.pred.ita$lower[7,]),
                          lmu = acct(qx.pred.ita$upper[7,]), lmx.s = qx.pred[1, ,7], lmi.s = qx.pred[2, ,7],
                          lmu.s = qx.pred[3, ,7], y = y[7,36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-5.15, -4.38)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt5

# 80+
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.ita$mean[13,]), lmi = acct(qx.m.ita$lower[13,]),
                     lmu = acct(qx.m.ita$upper[13,]), lmx.s = qx.sblc.rj$mean[13,], lmi.s = qx.sblc.rj$lower[13,],
                     lmu.s = qx.sblc.rj$upper[13,], y = y[13,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.ita$mean[13,]), lmi = acct(qx.pred.ita$lower[13,]),
                          lmu = acct(qx.pred.ita$upper[13,]), lmx.s = qx.pred[1, ,13], lmi.s = qx.pred[2, ,13],
                          lmu.s = qx.pred[3, ,13], y = y[13,36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-2.4, -1.85)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt6


## SMM
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.smm$mean[2,]), lmi = acct(qx.m.smm$lower[2,]),
                     lmu = acct(qx.m.smm$upper[2,]), lmx.s = qx.sblc.rj$mean[93,], lmi.s = qx.sblc.rj$lower[93,],
                     lmu.s = qx.sblc.rj$upper[93,], y = y[93,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.smm$mean[2,]), lmi = acct(qx.pred.smm$lower[2,]),
                          lmu = acct(qx.pred.smm$upper[2,]), lmx.s = qx.pred[1, ,93], lmi.s = qx.pred[2, ,93],
                          lmu.s = qx.pred[3, ,93], y = y[93,36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-6.6, -5.7)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt7


# 50-54
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.smm$mean[7,]), lmi = acct(qx.m.smm$lower[7,]),
                     lmu = acct(qx.m.smm$upper[7,]), lmx.s = qx.sblc.rj$mean[98,], lmi.s = qx.sblc.rj$lower[98,],
                     lmu.s = qx.sblc.rj$upper[98,], y = y[98,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.smm$mean[7,]), lmi = acct(qx.pred.smm$lower[7,]),
                          lmu = acct(qx.pred.smm$upper[7,]), lmx.s = qx.pred[1, ,98], lmi.s = qx.pred[2, ,98],
                          lmu.s = qx.pred[3, ,98], y = y[98,36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-5.15, -4.38)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt8

# 80+
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.smm$mean[13,]), lmi = acct(qx.m.smm$lower[13,]),
                     lmu = acct(qx.m.smm$upper[13,]), lmx.s = qx.sblc.rj$mean[104,], lmi.s = qx.sblc.rj$lower[104,],
                     lmu.s = qx.sblc.rj$upper[104,], y = y[104,1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.smm$mean[13,]), lmi = acct(qx.pred.smm$lower[13,]),
                          lmu = acct(qx.pred.smm$upper[13,]), lmx.s = qx.pred[1, ,104], lmi.s = qx.pred[2, ,104],
                          lmu.s = qx.pred[3, ,104], y = y[104,36:40])

ggplot() + theme_bw() + coord_cartesian(ylim = c(-2.65, -1.85)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt9


pdf("rj_rj_qxs.pdf", width=12, height=7)
pt1 + pt2 + pt3 + pt4 + pt5 + pt6 + pt7 + pt8 + pt9 + plot_layout(guides = "collect") &
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14))
graphics.off()

### cadeias



#### JAPAO ----
jp <- st_read("jpn_adm_2019_shp/jpn_admbnda_adm1_2019.shp")


aux.names <- c("Aichi", "Gifu", "Fukui", "Ishikawa", "Toyama", "Niigata",
                 "Nagano", "Yamanashi", "Shizuoka", "Gunma", "Tochigi", "Saitama",
                 "Ibaraki", "Chiba", "Tokyo", "Kanagawa", "Fukushima", "Yamagata")
aux.ages <- c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
              "55-59", "60-64", "65-69", "70-74", "75-79", "80-84")

jp %>%
  mutate(ADM1_EN = str_trim(ADM1_EN)) %>%
  filter(ADM1_EN %in% aux.names) -> jp2

pdf("jp_spatial.pdf", width=7, height=3)
ggplot(jp) + geom_sf() + geom_sf(data = jp2, aes(fill = ADM1_EN)) + theme_bw() + labs(fill = "") +
  coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) + guides(fill=guide_legend(ncol=2))
graphics.off()


### heatmap 25-29 em 1980-2000-2019
y <- matrix(NA, nrow = ages*n, ncol = t)
for(i in 1:n){
  aux <- aux.names[i]
  mx <- as.data.frame(read.table(paste0("japan_data/", str_trim(as.character(aux)), "_5x1.txt"), skip = 1, header = T))[,-c(3,4)] %>%
    filter(Year %in% 1980:2019, Age %in% aux.ages) %>%
    mutate(Total = log(as.numeric(Total))) %>%
    pivot_wider(names_from = Year, values_from = Total)
  y[((1:13) + 13*(i-1)),] <- as.matrix(mx[,-1])
}
ind <- rep(aux.names, each = 13)

age.it <- seq(2,234, by = 13)
data.aux <- cbind(jp2, x1 = y[age.it,1], x2 = y[age.it,21], x3 = y[age.it,40])
data.aux.rang <- cbind(t(y[age.it,1]), t(y[age.it,21]), t(y[age.it,40]))

library(ggplot2)
ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x1)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) + 
  ggtitle("1980") -> plt1
ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x2)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) +
  ggtitle("2000") -> plt2
ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x3)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) +
  ggtitle("2019") -> plt3

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

library(patchwork)
plt1 + plt2 + plt3 +
  plot_annotation(title = "Mortality rates for 25-29 age group, JP", theme = theme(plot.title = element_text(color="black", size=16))) +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        low = "yellow", high = "red") & 
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14)) -> pt
pdf("jp_heatmap_2529.pdf", width=12, height=3)
pt
graphics.off()

## 50-54
age.it <- seq(7,234, by = 13)
data.aux <- cbind(jp2, x1 = y[age.it,1], x2 = y[age.it,21], x3 = y[age.it,40])
data.aux.rang <- cbind(t(y[age.it,1]), t(y[age.it,21]), t(y[age.it,40]))

library(ggplot2)
ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x1)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) + 
  ggtitle("1980") -> plt4
ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x2)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) +
  ggtitle("2000") -> plt5
ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x3)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) +
  ggtitle("2019") -> plt6

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

library(patchwork)
plt4 + plt5 + plt6 +
  plot_annotation(title = "Mortality rates for 50-54 age group, JP", theme = theme(plot.title = element_text(color="black", size=16))) +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        high = "#00695C", low = "#8BC34A") & 
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14)) -> pt
pdf("jp_heatmap_5054.pdf", width=12, height=3)
pt
graphics.off()


## 80+
age.it <- seq(13,234, by = 13)
data.aux <- cbind(jp2, x1 = y[age.it,1], x2 = y[age.it,21], x3 = y[age.it,40])
data.aux.rang <- cbind(t(y[age.it,1]), t(y[age.it,21]), t(y[age.it,40]))

ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x1)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) + 
  ggtitle("1980") -> plt7
ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x2)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) +
  ggtitle("2000") -> plt8
ggplot(data = data.aux) + geom_sf(data = jp) + geom_sf(aes(fill = exp(as.numeric(x3)))) + labs(fill = "mx") +
  theme_bw() + coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) +
  ggtitle("2019") -> plt9

lims <- c(min(exp(data.aux.rang)), max(exp(data.aux.rang)))

plt7 + plt8 + plt9 +
  plot_annotation(title = "Mortality rates for 80+ age group, JP", theme = theme(plot.title = element_text(color="black", size=16))) +
  plot_layout(guides = "collect") & scale_fill_gradient(limits = lims,
                                                        low = "#26C6DA", high = "#9C27B0") & 
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14)) -> pt
pdf("jp_heatmap_80.pdf", width=12, height=3)
pt
graphics.off()

### parametros
fit.jp <- readRDS("jp_micro_fit_model1_full.RDS")
# alpha
alp <- t(apply(fit.jp$alpha.chain, 2, quantile, c(0.5, 0.025, 0.975)))
colnames(alp) <- c("fit", "lower", "upper")
ggplot(data = alp) + geom_line(aes(x = 1:13, y = fit)) + theme_bw() +
  geom_ribbon(aes(x = 1:13, ymax = upper, ymin = lower), fill = "grey40", alpha = .4) +
  labs(y = expression(alpha[x]), x = "Age group", fill = "") +
  scale_x_continuous(breaks = c(1, 5, 9, 13), labels = c("20-24", "40-44", "60-64" ,"80+")) -> pt1

# beta
alp <- t(apply(fit.jp$beta.chain, 2, quantile, c(0.5, 0.025, 0.975)))
colnames(alp) <- c("fit", "lower", "upper")
ggplot(data = alp) + geom_line(aes(x = 1:13, y = fit)) + theme_bw() +
  geom_ribbon(aes(x = 1:13, ymax = upper, ymin = lower), fill = "grey40", alpha = .4) +
  labs(y = expression(beta[x]), x = "Age group", fill = "") +
  scale_y_continuous(limits = c(0, 0.03)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = c(1, 5, 9, 13), labels = c("20-24", "40-44", "60-64" ,"80+")) -> pt2

# gamma
alp <- t(apply(fit.jp$gamma.chain, 2, quantile, c(0.5, 0.025, 0.975)))
colnames(alp) <- c("fit", "lower", "upper")
ggplot(data = alp) + geom_line(aes(x = 1:13, y = fit)) + theme_bw() +
  geom_ribbon(aes(x = 1:13, ymax = upper, ymin = lower), fill = "grey40", alpha = .4) +
  labs(y = expression(gamma[x]), x = "Age group", fill = "") +
  scale_x_continuous(breaks = c(1, 5, 9, 13), labels = c("20-24", "40-44", "60-64" ,"80+")) -> pt3

# kappa
alp <- t(apply(fit.jp$kappa.chain, 2, quantile, c(0.5, 0.025, 0.975)))
colnames(alp) <- c("fit", "lower", "upper")
ggplot(data = alp) + geom_line(aes(x = 1980:2014, y = fit)) + theme_bw() +
  geom_ribbon(aes(x = 1980:2014, ymax = upper, ymin = lower), fill = "grey40", alpha = .4) +
  labs(y = expression(kappa[t]), x = "Year", fill = "") -> pt4

# theta
alp <- apply(fit.rj$theta.chain, 2, median) ##menos diminui mort e positivo aumenta
data.aux <- cbind(jp2, theta = alp)
ggplot(data.aux) + geom_sf(aes(fill = theta)) + theme_bw() + labs(fill = "") +
  guides(fill = guide_colorbar(title = expression(theta[s]))) + 
  coord_sf(xlim = c(134,142.1), ylim = c(34.6,39.2)) +
  scale_fill_gradient2(low = "blue", midpoint = 0, high = "red") -> pt5

pdf("jp_parameters.pdf", width=12, height=4)
design <- "125
           345"
pt1 + pt2 + pt3 + pt4 + free(pt5, side = "tb") + plot_layout(guides = "collect", design = design) &
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14))
graphics.off()


### comp marginal
source("fitted_sblc.R")
qx.sblc.jp <- fitted_sblc(fit.jp)
source("predict_sblc.R")
#qx.sblc.pred <- predict.SBLC(fit.jp, h = 5)
#saveRDS(qx.sblc.pred, "jp_micro_qxpred.RDS")
qx.sblc.pred <- readRDS("jp_micro_qxpred.RDS")
qx.pred <- apply(qx.sblc.pred$y, 2:3, quantile, c(0.5, 0.025, 0.975))


library(BayesMortalityPlus)
## tokyo
tok <- which(ind == "Tokyo")
fit.mar.tk <- blc(y[tok, 1:35])
#plot(fit.mar.tk)
qx.m.tk <- fitted(fit.mar.tk)
qx.pred.tk <- fitted(predict(fit.mar.tk, h = 5))

## shizuoka
shi <- which(ind == "Shizuoka")
fit.mar.sh <- blc(y[shi, 1:35])
#plot(fit.mar.sh)
qx.m.sh <- fitted(fit.mar.sh)
qx.pred.sh <- fitted(predict(fit.mar.sh, h = 5))

# FUKUI
fk <- which(ind == "Fukui")
fit.mar.fk <- blc(y[fk, 1:35])
#plot(fit.mar.fk)
qx.m.fk <- fitted(fit.mar.fk)
qx.pred.fk <- fitted(predict(fit.mar.fk, h = 5))

acct <- function(x){log(-log(1 - x))}

## tokyo
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.tk$mean[2,]), lmi = acct(qx.m.tk$lower[2,]),
                     lmu = acct(qx.m.tk$upper[2,]), lmx.s = qx.sblc.jp$mean[tok[2],], lmi.s = qx.sblc.jp$lower[tok[2],],
                     lmu.s = qx.sblc.jp$upper[tok[2],], y = y[tok[2],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.tk$mean[2,]), lmi = acct(qx.pred.tk$lower[2,]),
                          lmu = acct(qx.pred.tk$upper[2,]), lmx.s = qx.pred[1, ,tok[2]], lmi.s = qx.pred[2, ,tok[2]],
                          lmu.s = qx.pred[3, ,tok[2]], y = y[tok[2],36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-6.6, -5.7)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "", title = "25-29") -> pt1


# 50-54
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.tk$mean[7,]), lmi = acct(qx.m.tk$lower[7,]),
                     lmu = acct(qx.m.tk$upper[7,]), lmx.s = qx.sblc.jp$mean[tok[7],], lmi.s = qx.sblc.jp$lower[tok[7],],
                     lmu.s = qx.sblc.jp$upper[tok[7],], y = y[tok[7],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.tk$mean[7,]), lmi = acct(qx.pred.tk$lower[7,]),
                          lmu = acct(qx.pred.tk$upper[7,]), lmx.s = qx.pred[1, ,tok[7]], lmi.s = qx.pred[2, ,tok[7]],
                          lmu.s = qx.pred[3, ,tok[7]], y = y[tok[7],36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-5.15, -4.38)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "", title = "50-54") -> pt2

# 80+
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.tk$mean[13,]), lmi = acct(qx.m.tk$lower[13,]),
                     lmu = acct(qx.m.tk$upper[13,]), lmx.s = qx.sblc.jp$mean[tok[13],], lmi.s = qx.sblc.jp$lower[tok[13],],
                     lmu.s = qx.sblc.jp$upper[tok[13],], y = y[tok[13],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.tk$mean[13,]), lmi = acct(qx.pred.tk$lower[13,]),
                          lmu = acct(qx.pred.tk$upper[13,]), lmx.s = qx.pred[1, ,tok[13]], lmi.s = qx.pred[2, ,tok[13]],
                          lmu.s = qx.pred[3, ,tok[13]], y = y[tok[13],36:40])

ggplot() + theme_bw() +# coord_cartesian(ylim = c(-2.4, -1.85)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "", title = "80+") -> pt3


### shizuoka
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.sh$mean[2,]), lmi = acct(qx.m.sh$lower[2,]),
                     lmu = acct(qx.m.sh$upper[2,]), lmx.s = qx.sblc.jp$mean[shi[2],], lmi.s = qx.sblc.jp$lower[shi[2],],
                     lmu.s = qx.sblc.jp$upper[shi[2],], y = y[shi[2],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.sh$mean[2,]), lmi = acct(qx.pred.sh$lower[2,]),
                          lmu = acct(qx.pred.sh$upper[2,]), lmx.s = qx.pred[1, ,shi[2]], lmi.s = qx.pred[2, ,shi[2]],
                          lmu.s = qx.pred[3, ,shi[2]], y = y[shi[2],36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-6.6, -5.7)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt4


# 50-54
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.sh$mean[7,]), lmi = acct(qx.m.sh$lower[7,]),
                     lmu = acct(qx.m.sh$upper[7,]), lmx.s = qx.sblc.jp$mean[shi[7],], lmi.s = qx.sblc.jp$lower[shi[7],],
                     lmu.s = qx.sblc.jp$upper[shi[7],], y = y[shi[7],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.sh$mean[7,]), lmi = acct(qx.pred.sh$lower[7,]),
                          lmu = acct(qx.pred.sh$upper[7,]), lmx.s = qx.pred[1, ,shi[7]], lmi.s = qx.pred[2, ,shi[7]],
                          lmu.s = qx.pred[3, ,shi[7]], y = y[shi[7],36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-5.15, -4.38)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt5

# 80+
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.sh$mean[13,]), lmi = acct(qx.m.sh$lower[13,]),
                     lmu = acct(qx.m.sh$upper[13,]), lmx.s = qx.sblc.jp$mean[shi[13],], lmi.s = qx.sblc.jp$lower[shi[13],],
                     lmu.s = qx.sblc.jp$upper[shi[13],], y = y[shi[13],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.sh$mean[13,]), lmi = acct(qx.pred.sh$lower[13,]),
                          lmu = acct(qx.pred.sh$upper[13,]), lmx.s = qx.pred[1, ,shi[13]], lmi.s = qx.pred[2, ,shi[13]],
                          lmu.s = qx.pred[3, ,shi[13]], y = y[shi[13],36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-2.4, -1.85)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt6


## Fukui
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.fk$mean[2,]), lmi = acct(qx.m.fk$lower[2,]),
                     lmu = acct(qx.m.fk$upper[2,]), lmx.s = qx.sblc.jp$mean[fk[2],], lmi.s = qx.sblc.jp$lower[fk[2],],
                     lmu.s = qx.sblc.jp$upper[fk[2],], y = y[fk[2],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.fk$mean[2,]), lmi = acct(qx.pred.fk$lower[2,]),
                          lmu = acct(qx.pred.fk$upper[2,]), lmx.s = qx.pred[1, ,fk[2]], lmi.s = qx.pred[2, ,fk[2]],
                          lmu.s = qx.pred[3, ,fk[2]], y = y[fk[2],36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-6.6, -5.7)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt7


# 50-54
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.fk$mean[7,]), lmi = acct(qx.m.fk$lower[7,]),
                     lmu = acct(qx.m.fk$upper[7,]), lmx.s = qx.sblc.jp$mean[fk[7],], lmi.s = qx.sblc.jp$lower[fk[7],],
                     lmu.s = qx.sblc.jp$upper[fk[7],], y = y[fk[7],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.fk$mean[7,]), lmi = acct(qx.pred.fk$lower[7,]),
                          lmu = acct(qx.pred.fk$upper[7,]), lmx.s = qx.pred[1, ,fk[7]], lmi.s = qx.pred[2, ,fk[7]],
                          lmu.s = qx.pred[3, ,fk[7]], y = y[fk[7],36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-5.15, -4.38)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt8

# 80+
aux.dt <- data.frame(year = 1980:2014, lmx = acct(qx.m.fk$mean[13,]), lmi = acct(qx.m.fk$lower[13,]),
                     lmu = acct(qx.m.fk$upper[13,]), lmx.s = qx.sblc.jp$mean[fk[13],], lmi.s = qx.sblc.jp$lower[fk[13],],
                     lmu.s = qx.sblc.jp$upper[fk[13],], y = y[fk[13],1:35])
aux.dt.pred <- data.frame(year = 2015:2019, lmx = acct(qx.pred.fk$mean[13,]), lmi = acct(qx.pred.fk$lower[13,]),
                          lmu = acct(qx.pred.fk$upper[13,]), lmx.s = qx.pred[1, ,fk[13]], lmi.s = qx.pred[2, ,fk[13]],
                          lmu.s = qx.pred[3, ,fk[13]], y = y[fk[13],36:40])

ggplot() + theme_bw() + #coord_cartesian(ylim = c(-2.65, -1.85)) +
  geom_point(data = aux.dt, aes(x = year, y = y, color = "Observed"), shape = 16) +
  geom_line(data = aux.dt, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu, ymin = lmi, color = "BLC"), linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt, aes(x = year, ymax = lmu.s, ymin = lmi.s, color = "SBLC"), linetype = "blank", fill = "red", alpha = .3) +
  
  geom_point(data = aux.dt.pred, aes(x = year, y = y, color = "Holdout"), shape = 4, stroke = 1) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx, color = "BLC")) +
  geom_line(data = aux.dt.pred, aes(x = year, y = lmx.s, color = "SBLC")) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu, ymin = lmi), color = "blue", linetype = "blank", fill = "blue", alpha = .3) +
  geom_ribbon(data = aux.dt.pred, aes(x = year, ymax = lmu.s, ymin = lmi.s), color = "red", linetype = "blank", fill = "red", alpha = .3) +
  scale_color_manual(values = c("SBLC" = "red", "BLC" = "blue", "Observed" = "grey4", "Holdout" = "grey4"),
                     breaks = c("Holdout", "Observed", "BLC", "SBLC")) +
  labs(y = "log(mx)", x = "Year", color = "", shape = "") -> pt9


pdf("jp_qxs.pdf", width=12, height=7)
pt1 + pt2 + pt3 + pt4 + pt5 + pt6 + pt7 + pt8 + pt9 + plot_layout(guides = "collect") &
  theme(legend.text=element_text(color="black", size=14),
        legend.title = element_text(color="black", size=14),
        axis.title = element_text(color = "black", size = 14))
graphics.off()
