## Writen by Lucas R. Moreira
## Usage: Calculates median values of LRT across sliding windows of 50 SNPs, whit 10 SNPs increments

library(data.table)
library(windowscanr)

PC1T <- fread("PC.res.nostruct.lrt0.gz")
PC2T <- fread("PC.res.nostruct.lrt1.gz")
PC3T <- fread("PC.res.nostruct.lrt2.gz")
PC1P <- fread("PC.res.nostruct.lrt3.gz")
PC2P <- fread("PC.res.nostruct.lrt4.gz")
PC3P <- fread("PC.res.nostruct.lrt5.gz")
PClist <- list(PC1T,PC2T,PC3T,PC1P,PC2P,PC3P)
names <- c("PC1T_nostruct","PC2T_nostruct","PC3T_nostruct","PC1P_nostruct","PC2P_nostruct","PC3P_nostruct")

for(i in 1:length(PClist)){
  hs <- PClist[[i]]
  name <- names[i]
  x <- hs[hs$beta!="NaN",]
  rol_win <- winScan(x = x,
                     groups = "Chromosome",
                     position = NULL,
                     values = "LRT",
                     win_size = 50,
                     win_step = 10,
                     funs = "median",
                     cores = 32)
  write.table(rol_win,paste0(name,".rol_win"),row.names = F)
}
