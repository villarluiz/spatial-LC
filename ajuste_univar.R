library(tidyverse)
library(BayesMortalityPlus)

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
plot(fit.smm, parameter = "fitted", age = as.integer(seq(20,80, by = 5)))
## satisfatório

#### serie temporal mortalidade x anos
qxs <- fitted(fit)
aux.qx <- matrix(nrow = 3, ncol = 42)
aux.qx[1,] <- apply(qxs$lower, 2, mean)
aux.qx[2,] <- apply(qxs$mean, 2, mean)
aux.qx[3,] <- apply(qxs$upper, 2, mean)

plot(1980:2021, aux.qx[2,], lwd=2, type="l", ylim = c(min(aux.qx), max(aux.qx)))
lines(1980:2021, aux.qx[1,], lwd = 2, lty=2, col = "blue")
lines(1980:2021, aux.qx[3,], lwd = 2, lty=2, col = "blue")

plot(fit, parameter = "kappa")


# modelar vizinhos de sta maria madalena para comparacao ----
# Santo antonio de padua ## 33002 
y <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt"))
  y[,i] <- c(as.numeric(mx$X33002.SANTO.ANTONIO.DE.PADUA))
}
colnames(y) <- as.character(1980:2021)
row.names(y) <- mx[,1]


set.seed(14)
fit.sap <- blc(y, M = 10000)
plot(fit.sap)
plot(fit.sap, parameter = "fitted", age = as.integer(seq(20,80, by = 5)))
## satisfatório


# Campo dos goytacazes  ## 33003 
y <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt"))
  y[,i] <- c(as.numeric(mx$X33003.CAMPOS.DOS.GOYTACAZES))
}
colnames(y) <- as.character(1980:2021)
row.names(y) <- mx[,1]

set.seed(14)
fit.cg <- blc(y, M = 10000)
plot(fit.cg)
plot(fit.cg, parameter = "fitted", age = as.integer(seq(20,80, by = 5)))
## satisfatório

# Macae ## 33004
y <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt"))
  y[,i] <- c(as.numeric(mx$X33004.MACAE))
}
colnames(y) <- as.character(1980:2021)
row.names(y) <- mx[,1]


set.seed(14)
fit.mac <- blc(y, M = 10000)
plot(fit.mac)
plot(fit.mac, parameter = "fitted", age = as.integer(seq(20,80, by = 5)))
## satisfatório

# Cantagalo-Cordeiro ## 33006
y <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt"))
  y[,i] <- c(as.numeric(mx$X33006.CANTAGALO.CORDEIRO))
}
colnames(y) <- as.character(1980:2021)
row.names(y) <- mx[,1]


set.seed(14)
fit.cc <- blc(y, M = 10000)
plot(fit.cc)
plot(fit.cc, parameter = "fitted", age = as.integer(seq(20,80, by = 5)))
## satisfatório

# Nova Friburgo ## 33007 
y <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt"))
  y[,i] <- c(as.numeric(mx$X33007.NOVA.FRIBURGO))
}
colnames(y) <- as.character(1980:2021)
row.names(y) <- mx[,1]


set.seed(14)
fit.nv <- blc(y, M = 10000)
plot(fit.nv)
plot(fit.nv, parameter = "fitted", age = as.integer(seq(20,80, by = 5)))
## satisfatório

library(patchwork)
library(ggplot2)
### beta
aux = "beta"
plot(fit.smm, parameter = aux) + ggtitle("Sta. Maria Madalena") +
  plot(fit.cc, parameter = aux) + ggtitle("Cantagalo-Cordeiro") +
  plot(fit.cg, parameter = aux) + ggtitle("Campo dos Goytacazes") +
  plot(fit.mac, parameter = aux) + ggtitle("Macae") +
  plot(fit.nv, parameter = aux) + ggtitle("Nova Friburgo") +
  plot(fit.sap, parameter = aux) + ggtitle("St. Antonio de Padua") -> pt.b
ggsave("graficos/comp_betas.png", pt.b, device = png(width = 600, height = 400))
graphics.off()

### kappa
aux = "kappa"
plot(fit.smm, parameter = aux) + ggtitle("Sta. Maria Madalena") +
  plot(fit.cc, parameter = aux) + ggtitle("Cantagalo-Cordeiro") +
  plot(fit.cg, parameter = aux) + ggtitle("Campo dos Goytacazes") +
  plot(fit.mac, parameter = aux) + ggtitle("Macae") +
  plot(fit.nv, parameter = aux) + ggtitle("Nova Friburgo") +
  plot(fit.sap, parameter = aux) + ggtitle("St. Antonio de Padua") -> pt.k
ggsave("graficos/comp_kappas.png", pt.k, device = png(width = 600, height = 400))
graphics.off()

### comparar alphas no mesmo grafico
x = as.integer(seq(20,80, by = 5))
ggplot(NULL) +
  geom_line(data = data.frame(y = apply(fit.smm$alpha, 1, mean), x = x), aes(x = x, y = y, col = "Sta. Maria Madalena")) +
  geom_ribbon(data = data.frame(t(rbind(apply(fit.smm$alpha, 1, quantile, c(0.025, 0.975)), x))),
              aes(x = x, ymin = X2.5., ymax = X97.5., fill = "Sta. Maria Madalena IC"), alpha = 0.2) +
  
  geom_line(data = data.frame(y = apply(fit.cc$alpha, 1, mean), x = x), aes(x = x, y = y, col = "Cantagalo-Cordeiro")) +
  geom_ribbon(data = data.frame(t(rbind(apply(fit.cc$alpha, 1, quantile, c(0.025, 0.975)), x))),
              aes(x = x, ymin = X2.5., ymax = X97.5., fill = "Cantagalo-Cordeiro IC"), alpha = 0.2) +
  
  geom_line(data = data.frame(y = apply(fit.cg$alpha, 1, mean), x = x), aes(x = x, y = y, col = "Campo dos Goytacazes")) +
  geom_ribbon(data = data.frame(t(rbind(apply(fit.cg$alpha, 1, quantile, c(0.025, 0.975)), x))),
              aes(x = x, ymin = X2.5., ymax = X97.5., fill = "Campo dos Goytacazes IC"), alpha = 0.2) +
  
  geom_line(data = data.frame(y = apply(fit.mac$alpha, 1, mean), x = x), aes(x = x, y = y, col = "Macae")) +
  geom_ribbon(data = data.frame(t(rbind(apply(fit.mac$alpha, 1, quantile, c(0.025, 0.975)), x))),
              aes(x = x, ymin = X2.5., ymax = X97.5., fill = "Macae IC"), alpha = 0.2) +
  
  geom_line(data = data.frame(y = apply(fit.nv$alpha, 1, mean), x = x), aes(x = x, y = y, col = "Nova Friburgo")) +
  geom_ribbon(data = data.frame(t(rbind(apply(fit.nv$alpha, 1, quantile, c(0.025, 0.975)), x))),
              aes(x = x, ymin = X2.5., ymax = X97.5., fill = "Nova Friburgo IC"), alpha = 0.2) +
  
  geom_line(data = data.frame(y = apply(fit.sap$alpha, 1, mean), x = x), aes(x = x, y = y, col = "St. Antonio de Padua")) +
  geom_ribbon(data = data.frame(t(rbind(apply(fit.sap$alpha, 1, quantile, c(0.025, 0.975)), x))),
              aes(x = x, ymin = X2.5., ymax = X97.5., fill = "St. Antonio de Padua IC"), alpha = 0.2) +
  
  xlab("age group") + ylab("alpha") + theme_bw() +
  scale_color_manual(values = rainbow(6)) +
  scale_fill_manual(values = rainbow(6)) +
  theme(axis.title.x = ggplot2::element_text(color = 'black', size = 13),
        axis.title.y = ggplot2::element_text(color = 'black', size = 13)) -> g
plotly:::ggplotly(g)

# modelagem univariado, RJ (bom volume de dados) ##33018 -----
y <- matrix(NA, nrow = 13, ncol = 42)
for(i in 1:42){
  aux <- i + 1979
  mx <- read.table(paste0("Dados_microrregioes/det/logmx_det/logmx_", as.character(aux), ".txt"))
  y[,i] <- c(as.numeric(mx$X33018.RIO.DE.JANEIRO))
}
colnames(y) <- as.character(1980:2021)
row.names(y) <- mx[,1]


set.seed(14)
fit <- blc(y[4:13,], M = 10000)
plot(fit)
plot(fit, parameter = "fitted", age = as.integer(seq(35,80, by = 5)))
