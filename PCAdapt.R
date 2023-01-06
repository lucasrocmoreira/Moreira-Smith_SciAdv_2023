#############################################
#                  PCAdapt                  #
#              Lucas R. Moreira             #
#############################################

#install.packages("pcadapt")
#devtools::install_github('kaustubhad/fastman',build_vignettes = TRUE)
library(RcppCNPy)
library(bigutilsr)
library(fastman)
library(data.table)

zscores <- npyLoad("PP.pcangsd.pcadapt.zscores.npy")
K <- ncol(zscores)

sites <- fread("Picoides-pubescens.positions")
filtered.sites <- fread("PP.pcangsd.sites")
sites <- sites[filtered.sites$V1==1,]
sites$snp_id <- paste0(sites$chr,"_",sites$pos)

# Which are the sites also present in ANGSDasso?
ANGSDasso <- fread("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\Window-based\\sites.tsv")
ANGSDasso$snp_id <- paste0(ANGSDasso$Chromosome,"_",ANGSDasso$Position)

# For one component only
if (K == 1) {
  d2 <- (zscores - median(zscores))^2
} else {
  d2 <- dist_ogk(zscores)
}

p.values <- pchisq(d2, df=K, lower.tail=F)
subselect.p.value <- p.values[sites$snp_id %in% ANGSDasso$snp_id]
subselect.snp.position <- sites[sites$snp_id %in% ANGSDasso$snp_id,]
final.file.to.export <- cbind(subselect.snp.position,p.value=subselect.p.value)

write.table(final.file.to.export,"PP.p.values.tsv",quote=F, row.names=F)

write.table(d2, file="PP.pcadapt.test.txt", quote=F, row.names=F, col.names=F)
write.table(p.values, file="PP.pcadapt.pval.txt", quote=F, row.names=F, col.names=F)

################ Visualizing the data ################

hist(p.values, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

fastqq(p.values, speedup=TRUE, lambda=TRUE, fix_zero=TRUE, cex=0.6, cex.axis=0.9)

################ Choosing a cutoff for outlier detection ################

PP.roll <- fread("PP.rol_win")

# autosomal
quantile(PP.roll[PP.roll$chr!="Z",]$p.value_median,1-.995)
significance.auto.p <- 0.01

# Z
quantile(PP.roll[PP.roll$chr=="Z",]$p.value_median,1-.995)
significance.Z.p <- 1e-05

signficant.auto <- PP.roll[PP.roll$chr!="Z" & PP.roll$p.value_median<significance.auto.p,]
signficant.Z <- PP.roll[PP.roll$chr=="Z" & PP.roll$p.value_median<significance.Z.p,]

signficant <- rbind(signficant.auto,signficant.Z)

saveRDS(signficant,"PP.candidates.RDS")
write.table(signficant,"PP.candidates.tsv",row.names = F, quote = F, sep = '\t')

# get candidates windows for later comparison
candidates <- readRDS("PP.candidates.RDS")

write.csv(paste0(candidates$chr,"_",candidates$win_mid),"all.candidates.windows.csv",row.names = F,quote =F)
