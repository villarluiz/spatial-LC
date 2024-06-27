library(tidyverse)

### TESTE
ex <- as.data.frame(t(read.table("ex_det/ex_1996.txt", sep = ";")))
dx <- as.data.frame(t(read.table("obt_det/obt_1996_alt.txt", sep = ";")))
age <- dx[,1]; reg <- dx[1,]
aux.dx <- as.data.frame(sapply(dx[-1,2:64], as.numeric))
aux.ex <- as.data.frame(sapply(ex[-1,2:64], as.numeric))
logmx <- log(aux.dx / aux.ex)[-14,]

plot(logmx[,1], type = "l", ylim = c(-10,-1))
for(i in 2:18){
  lines(logmx[,i], col = "red")
}  ## OK

logmx <- data.frame(age = cbind(age[-c(1,15)], logmx)) 
colnames(logmx) <- c("age", reg[1,-c(1,65)])
###salvar
write.table(logmx, "logmx_det/logmx_1996_det.txt")              


### FUNCAO
nas <- rep(NA, length(1980:2021))
for(i in seq(1980,2021)){
  aux.na <- 0
  ex <- as.data.frame(t(read.table(paste("ex_det/ex_", as.character(i), ".txt", sep = ""),sep = ";")))
  dx <- as.data.frame(t(read.table(paste("obt_det/obt_", as.character(i), ".txt", sep = ""),sep = ";")))
  age <- dx[,1]; reg <- dx[1,]
  aux.dx <- as.data.frame(sapply(dx[-1,2:64], as.numeric))
  aux.ex <- as.data.frame(sapply(ex[-1,2:64], as.numeric))
  logmx <- log(aux.dx / aux.ex)[-14,]
  nas[i-1979] <- sum(is.na(logmx))
  
  logmx <- data.frame(age = cbind(age[-c(1,15)], logmx)) 
  colnames(logmx) <- c("age", reg[1,-c(1,65)])
  write.table(logmx, paste("logmx_det/logmx_", as.character(i), ".txt", sep = ""))
}

### analise do intervalo de idades a ser usado puro
rbind(1980:2021, nas) ##1985, 1991, 1992, 2002, 2018, 2021
data <- read.table("logmx_det/logmx_1985.txt")
data <- read.table("logmx_det/logmx_2002.txt")
data <- read.table("logmx_det/logmx_2018.txt")
data <- read.table("logmx_det/logmx_2021.txt")
