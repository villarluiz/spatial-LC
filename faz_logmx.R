library(tidyverse)

### TESTE
ex <- as.data.frame(t(read.table("Dados_microrregioes/ex_OPS/ex_2021.txt")))
dx <- as.data.frame(read.table("Dados_microrregioes/obt_OPS/obt_2021.txt", sep = ";"))[,-c(12,13)]
reg <- dx[,1]
aux.dx <- as.data.frame(sapply(dx[-1,2:11], as.numeric))
logmx <- t(log(aux.dx / ex[-1,])) 

plot(logmx[,1], type = "l", ylim = c(-10,-1))
for(i in 2:19){
  lines(logmx[,i], col = "red")
}  ## OK

logmx <- data.frame(rbind(reg[-1], logmx)) 
###salvar
write.table(logmx, "Dados_microrregioes/logmx_OPS/logmx_2021_OPS.txt")              


### FUNCAO
for(i in seq(2021,2011)){
  ex <- as.data.frame(t(read.table(paste("Dados_microrregioes/ex_OPS/ex_", as.character(i), ".txt", sep = ""))))
  dx <- as.data.frame(read.table(paste("Dados_microrregioes/obt_OPS/obt_", as.character(i), ".txt", sep = ""),sep = ";"))[,-c(12,13)]
  reg <- dx[,1]
  aux.dx <- as.data.frame(sapply(dx[-1,2:11], as.numeric))
  logmx <- t(log(aux.dx / ex[-1,])) 
  
  logmx <- data.frame(rbind(reg[-1], logmx)) 
  write.table(logmx, paste("Dados_microrregioes/logmx_OPS/logmx_", as.character(i), ".txt", sep = ""))
}
