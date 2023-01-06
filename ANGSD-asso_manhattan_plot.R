## Writen by Lucas R. Moreira
## Usage: Plots sliding windows Manhattan plots

library(data.table)
library(fastman)
library(windowscanr)

# Read results
PC1T <- fread("PC1T.rol_win")  #function to read very large files
PC2T <- fread("PC2T.rol_win")
PC3T <- fread("PC3T.rol_win")
PC1P <- fread("PC1P.rol_win")
PC2P <- fread("PC2P.rol_win")
PC3P <- fread("PC3P.rol_win")
PClist <- list(PC1T,PC2T,PC3T,PC1P,PC2P,PC3P)
names <- c("PC1T","PC2T","PC3T","PC1P","PC2P","PC3P")

for(i in 1:length(PClist)){
  hs <- PClist[[i]]
  name <- names[i]
  print(name)
  print("Reordering and replacing the chromosome names by numbers")
  chrOrder<-c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","LG2","LGE22","Z","w")
  hs$Chromosome<-factor(hs$Chromosome, levels=chrOrder)
  df <- hs[order(hs$Chromosome),]
  
  x <- df
  x$Chromosome <- as.character(x$Chromosome)
  pdf(paste0(name,".window.LRT.pdf"),width = 10,height = 5)
  fastman(x, chr = "Chromosome", bp = "win_mid", p = "LRT_median",logp = F, sortchr = F, maxP = max(x$LRT_median), 
          speedup=T, col="rgbs",cex=0.5,xlab="Pseudochromosome",ylab="Median LRT",cex.axis = 1)
  dev.off()
}

# Plot PC2 and PC3 for Supplementary Information

PC.2.3.list <- list(PC2T,PC3T,PC2P,PC3P)
names.2.3 <- c("PC2T","PC3T","PC2P","PC3P")

pdf("PC2&3.window.LRT.pdf",width = 15,height = 10)
par(mfrow = c(2,2))
for(i in 1:length(PC.2.3.list)){
  hs <- PC.2.3.list[[i]]
  name <- names.2.3[i]
  print(name)
  print("Reordering and replacing the chromosome names by numbers")
  chrOrder<-c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","LG2","LGE22","Z","w")
  hs$Chromosome<-factor(hs$Chromosome, levels=chrOrder)
  df <- hs[order(hs$Chromosome),]
  
  x <- df
  x$Chromosome <- as.character(x$Chromosome)
  fastman(x, chr = "Chromosome", bp = "win_mid", p = "LRT_median",logp = F, sortchr = F, maxP = max(x$LRT_median), 
          speedup=T, col="rgbs",cex=0.5,xlab="Pseudochromosome",ylab="Median LRT",cex.axis = 1,main=paste("Downy -",name))
}
dev.off()