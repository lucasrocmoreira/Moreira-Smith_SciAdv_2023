## Writen by Lucas R. Moreira
## Usage: Genotype-environment association analysis

library(data.table)
library(dplyr)
library(ggplot2)
library(fastman)

ann <- fread("../PP-SNP_annotation.csv",header=T,sep = ",")
sites <- fread("sites.tsv")
which.sites <- paste0(sites$Chromosome,"_",sites$Position)
ann.sites <- ann[ann$snp_id %in% which.sites,]

# read the results from ANGSD-asso, considering the median values across sliding windows
PC1T <- fread("PC1T.rol_win")  #function to read very large files
PC2T <- fread("PC2T.rol_win")
PC3T <- fread("PC3T.rol_win")
PC1P <- fread("PC1P.rol_win")
PC2P <- fread("PC2P.rol_win")
PC3P <- fread("PC3P.rol_win")

# read the significance thresholds, considering the median values across sliding windows
LRT_PC1T_threshold_0.00001 <- fread("PC1T.ran.rol_win",header = T)
LRT_PC2T_threshold_0.00001 <- fread("PC2T.ran.rol_win",header = T)
LRT_PC3T_threshold_0.00001 <- fread("PC3T.ran.rol_win",header = T)
LRT_PC1P_threshold_0.00001 <- fread("PC1P.ran.rol_win",header = T)
LRT_PC2P_threshold_0.00001 <- fread("PC2P.ran.rol_win",header = T)
LRT_PC3P_threshold_0.00001 <- fread("PC3P.ran.rol_win",header = T)

file_list <- list(PC1T,PC2T,PC3T,PC1P,PC2P,PC3P)
names(file_list) <- c("PC1T","PC2T","PC3T","PC1P","PC2P","PC3P")
threshold_list <- list(LRT_PC1T_threshold_0.00001,LRT_PC2T_threshold_0.00001,
                       LRT_PC3T_threshold_0.00001,LRT_PC1P_threshold_0.00001,
                       LRT_PC2P_threshold_0.00001,LRT_PC3P_threshold_0.00001)

#Check distribution  of p-values:

find.p.value <- function(x){
  pchisq(x, df=1, lower.tail=FALSE)
}
PC1T_p.value <- unlist(lapply(PC1T$LRT_median,FUN = find.p.value))

pdf("qqplot.pdf", width=6, height=6)
fastqq(PC1T_p.value, speedup=TRUE, lambda=TRUE, fix_zero=TRUE, cex=0.6, cex.axis=0.9)
dev.off()

# The QQ plot is a graphical representation of the deviation of the observed P values from the null hypothesis: 
# the observed P values for each SNP are sorted from largest to smallest and plotted against expected values from a theoretical χ2-distribution. 
# If the observed values correspond to the expected values, all points are on or near the middle line between the x-axis and the y-axis.
# If some observed P values are clearly more significant than expected under the null hypothesis, points will move towards the y-axis.
# If there is an early separation of the expected from the observed, this means that many moderately significant P values are more significant than expected under the null hypothesis.
# Reasons:
# Population structure
# low N

for(i in 1:length(file_list)){
  SNPs <- file_list[[i]]
  threshold <- threshold_list[[i]]
  name <- names(file_list)[i]
  print(name)
  
  # Alternative thresholds
  # 1)
  # 0.99999 from null distribution / error rate = 1e-06
  # FDR at 0.001 = 28
  #significant <- SNPs[SNPs$LRT>20,]
  
  # 2)
  # Top 0.001% of FDRs
  # quantile(SNPs$LRT,0.999)
  # significant <- SNPs[SNPs$LRT>14,]
  
  # 3) 
  # Above 0.0001 of random permutation
  significant <- SNPs[SNPs$LRT_median>threshold$LRT_median,]
  print(paste("Saving",name))
  saveRDS(significant,file=paste0(name,".RDS"))
  write.csv(significant,paste0("Candidates_info_",name,".csv"),row.names = F,quote = F)
}


######################### INVESTIGATE ANNOTATION ############################

candidates_PC1T <- readRDS("PC1T.RDS")
candidates_PC2T <- readRDS("PC2T.RDS")
candidates_PC3T <- readRDS("PC3T.RDS")
candidates_PC1P <- readRDS("PC1P.RDS")
candidates_PC2P <- readRDS("PC2P.RDS")
candidates_PC3P <- readRDS("PC3P.RDS")
env.variables <- list(candidates_PC1T,candidates_PC2T,candidates_PC3T,
                      candidates_PC1P,candidates_PC2P,candidates_PC3P)
names(env.variables) <- c("PC1T","PC2T","PC3T","PC1P","PC2P","PC3P")

# get candidates for temp and precip
candidates_PC1T$window_id <- paste0(candidates_PC1T$Chromosome,"_",candidates_PC1T$win_mid)
candidates_PC2T$window_id <- paste0(candidates_PC2T$Chromosome,"_",candidates_PC2T$win_mid)
candidates_PC3T$window_id <- paste0(candidates_PC3T$Chromosome,"_",candidates_PC3T$win_mid)
candidates_PC1P$window_id <- paste0(candidates_PC1P$Chromosome,"_",candidates_PC1P$win_mid)
candidates_PC2P$window_id <- paste0(candidates_PC2P$Chromosome,"_",candidates_PC2P$win_mid)
candidates_PC3P$window_id <- paste0(candidates_PC3P$Chromosome,"_",candidates_PC3P$win_mid)

all_candidates <- unique(c(candidates_PC1T$window_id,candidates_PC2T$window_id,candidates_PC3T$window_id,
                           candidates_PC1P$window_id,candidates_PC2P$window_id,candidates_PC3P$window_id))
write.csv(all_candidates,"all.candidates.windows.csv",row.names = F,quote =F)

# Loop through candidates for different variables
for(i in 1:length(env.variables)){
  candidates_PC <- env.variables[[i]]
  name <- names(env.variables)[i]
  print(name)
  
  candidate.genes <- numeric()
  annotation.genes <- data.frame()
  for(i in 1:nrow(candidates_PC)){
    print(paste0("Window ",i," of ",nrow(candidates_PC)))
    chrom <- candidates_PC[i,]$Chromosome
    start <- candidates_PC[i,]$win_start
    end <- candidates_PC[i,]$win_end
    LRT <- candidates_PC[i,]$LRT_median
    ann.chr <- ann.sites[ann.sites$chromosome.x==chrom,]
    annotation <- ann.chr[start:end,]
    annotation$LRT <- LRT
    annotation.genes <- rbind(annotation.genes,annotation)
    annotation.genes <- annotation.genes %>% distinct()
    
    genes_ann <- annotation$ID
    genes <- as.data.frame(numeric())
    for(i in 1:length(genes_ann)){
      snp.gene <- genes_ann[i]
      if(grepl("-",snp.gene,fixed=T)){
        split.gene <- unlist(strsplit(snp.gene,"-"))
        genes <- unique(c(genes,split.gene))
      }else{
        genes <- unique(c(genes,snp.gene))
      }
    }
    candidate.genes <- unlist(unique(c(candidate.genes,genes)))
  }
  
  write.csv(candidate.genes,paste0(name,"_gene_id.csv"),row.names = F)
  write.csv(annotation.genes,paste0(name,"_snp_annotation.csv"),row.names = F)
  
}  

# How many unique windows for temp and precip?
uniq.temp <- unique(c(paste0(candidates_PC1T$Chromosome,"_",candidates_PC1T$win_mid),
                    paste0(candidates_PC2T$Chromosome,"_",candidates_PC2T$win_mid),
                    paste0(candidates_PC3T$Chromosome,"_",candidates_PC3T$win_mid)))

uniq.prec <- unique(c(paste0(candidates_PC1P$Chromosome,"_",candidates_PC1P$win_mid),
                      paste0(candidates_PC2P$Chromosome,"_",candidates_PC2P$win_mid),
                      paste0(candidates_PC3P$Chromosome,"_",candidates_PC3P$win_mid)))

  
################ GO TERM ENRICHMENT ################

candidates_PC1T <- read.csv("PC1T_gene_id.csv")$x
candidates_PC2T <- read.csv("PC2T_gene_id.csv")$x
candidates_PC3T <- read.csv("PC3T_gene_id.csv")$x
candidates_PC1P <- read.csv("PC1P_gene_id.csv")$x
candidates_PC2P <- read.csv("PC2P_gene_id.csv")$x
candidates_PC3P <- read.csv("PC3P_gene_id.csv")$x
candidates_T <- c(candidates_PC1T,candidates_PC2T,candidates_PC3T)
candidates_P <- c(candidates_PC1P,candidates_PC2P,candidates_PC3P)

temp.all.gene_ids <- unique(candidates_T)
prec.all.gene_ids <- unique(candidates_P)

write.csv(temp.all.gene_ids,"temp.all.gene_ids.csv",row.names = F)
write.csv(prec.all.gene_ids,"prec.all.gene_ids.csv",row.names = F)

go <- read.delim("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\PANNZER2\\PANNZER2 Gene Ontology\\GO.out.txt",sep="\t",header=TRUE)

####################################
#                                 
#
#            TEMPERATURE
#
#
###################################

# Here I'm considering only genes containing SNPs within or nearby
temp.genes <- read.csv("temp.all.gene_ids.csv")$x

total_number_of_genes_in_genome <- length(unique(go$qpid))
total_number_of_candidates <- length(temp.genes)
total_genes <- cbind(total_number_of_candidates,total_number_of_genes_in_genome)

BP <- numeric()
for(i in 1:length(temp.genes)){
  candidate <- temp.genes[i]
  print(candidate)
  go_set <- go[grep(candidate,go$qpid),]
  #go_set
  go_BP <- go_set[grep("BP",go_set$ontology),]
  BP <- rbind(BP,go_BP)
}
#write.csv(BP,"GO.TEMP.BP.only_within.csv",row.names = F)

###################################
# Enrichment analysis (MANUALLY)  #
###################################

# BIOLOGICAL PROCESS

# Count how many genes in the candidate set have each GO id
candidates_goid <- as.factor(sprintf("%07d", BP$goid))
t.candidates_goid <- table(candidates_goid)

# Count how many genes in the entire annotation have each GO id
go.genome <- go[grep("BP",go$ontology),]
genome_goid <- as.factor(sprintf("%07d", go.genome$goid))
matches <- droplevels(genome_goid[is.element(genome_goid, candidates_goid)])
t.genome_goid <- table(matches)

# Produce a contingency table
contingency <- cbind(candidate_set=t.candidates_goid,genome=(t.genome_goid-t.candidates_goid))

######################
# Fisher exact test  #
######################

fisher_test <- numeric()
for(i in 1:nrow(contingency)){
  target <- contingency[i,]
  not.target <- total_genes-target
  contingency.table <- rbind(target=target,others=not.target)
  result <- fisher.test(contingency.table,alternative = "greater")
  estimate <- unname(result$estimate)
  p <- result$p.value
  go_id <- paste0("GO:",rownames(contingency)[i])
  row.info <- cbind("GO_ID"=go_id,"candidate_count"=unname(target[1]),"genome_count"=unname(target[1]+target[2]),"odds.ratio"=estimate,"p.value"=p)
  fisher_test <- rbind(fisher_test,row.info)
}
fisher_test <- fisher_test[order(fisher_test[,5]),]

corrected.p <- p.adjust(fisher_test[,5], method = "fdr")
fisher_test <- as.data.frame(cbind(fisher_test,corrected.p))
write.csv(fisher_test,"BP-fisher_test.enrichment.csv",row.names = F)

enrichment <- fisher_test[which(as.numeric(fisher_test$corrected.p)<0.05),]

write.csv(enrichment,"BP-temp-enrichment.csv",row.names = F)

####################################
#                                 
#
#            PRECIPITATION
#
#
###################################

# Here I'm considering only genes containing SNPs within or nearby
prec.genes <- read.csv("prec.all.gene_ids.csv")$x

total_number_of_genes_in_genome <- length(unique(go$qpid))
total_number_of_candidates <- length(prec.genes)
total_genes <- cbind(total_number_of_candidates,total_number_of_genes_in_genome)

BP <- numeric()
for(i in 1:length(prec.genes)){
  candidate <- prec.genes[i]
  print(candidate)
  go_set <- go[grep(candidate,go$qpid),]
  #go_set
  go_BP <- go_set[grep("BP",go_set$ontology),]
  BP <- rbind(BP,go_BP)
}
#write.csv(BP,"GO.PREC.BP.only_within.csv",row.names = F)

###################################
# Enrichment analysis (MANUALLY)  #
###################################

# BIOLOGICAL PROCESS

# Count how many genes in the candidate set have each GO id
candidates_goid <- as.factor(sprintf("%07d", BP$goid))
t.candidates_goid <- table(candidates_goid)

# Count how many genes in the entire annotation have each GO id
go.genome <- go[grep("BP",go$ontology),]
genome_goid <- as.factor(sprintf("%07d", go.genome$goid))
matches <- droplevels(genome_goid[is.element(genome_goid, candidates_goid)])
t.genome_goid <- table(matches)

# Produce a contingency table
contingency <- cbind(candidate_set=t.candidates_goid,genome=(t.genome_goid-t.candidates_goid))

######################
# Fisher exact test  #
######################

fisher_test <- numeric()
for(i in 1:nrow(contingency)){
  target <- contingency[i,]
  not.target <- total_genes-target
  contingency.table <- rbind(target=target,others=not.target)
  result <- fisher.test(contingency.table)
  estimate <- unname(result$estimate)
  p <- result$p.value
  go_id <- paste0("GO:",rownames(contingency)[i])
  row.info <- cbind("GO_ID"=go_id,"candidate_count"=unname(target[1]),"genome_count"=unname(target[1]+target[2]),"odds.ratio"=estimate,"p.value"=p)
  fisher_test <- rbind(fisher_test,row.info)
}

corrected.p <- p.adjust(fisher_test[,5], method = "fdr")
fisher_test <- as.data.frame(cbind(fisher_test,corrected.p))
write.csv(fisher_test,"BP-prec-fisher_test.enrichment.csv",row.names = F)

enrichment <- fisher_test[which(as.numeric(fisher_test$corrected.p)<0.05),]

write.csv(enrichment,"BP-prec-enrichment.csv",row.names = F)
