## Made by Lucas R. Moreira
## Last updated 28 March 2020
## Usage: Calculates genotype-environemnt association

# install.packages("devtools")
# devtools::install_github("bcm-uga/lfmm")

#BiocManager::install("LEA")

library(lfmm)
library(data.table)
library(dplyr)
library(ggplot2)

# read genotypes

lfmm <- fread("PP.SNP-only_simplified.imputed.lfmm")  #function to read very large files
Y <- lfmm

# ################
# # running snmf #
# ################
# 
# project.snmf = snmf("SNP-only_simplified.lfmm", K = 4, 
#                     entropy = TRUE, repetitions = 10,
#                     project = "new")
# 
# # select the run with the lowest cross-entropy value
# best = which.min(cross.entropy(project.snmf, K = 4))
# 
# # Impute the missing genotypes
# impute(project.snmf, "SNP-only_simplified.lfmm", method = 'mode', K = 4, run = best)


# read environmental data

env <- read.csv("Picoides_pubescens_PCscores.csv",header=T)
env_data <- env[2:7]
X <- env_data

# run model

mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 4)  # K = 4 genetic clusters

pp <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif")

pvalues <- pp$calibrated.pvalue

#saveRDS(pvalues,"pvalues_Temp.RDS")
pvalues <- readRDS("pvalues_Temp.RDS")

ann <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\PP-SNP_annotation.csv",header=T)

####################################
# Find candidate SNPs
####################################

for(i in 1:ncol(pvalues)){
  PC <- pvalues[,i]
  name <- colnames(pvalues)[i]
  print(name)

  x <- as.data.frame(PC)
  
  PC_reduced <- x %>%
    filter(-log10(PC)>1.5)  # Reduce the number of SNPs for plotting
  
  ## FDR control: Benjamini-Hochberg at level q
  ## L = number of loci
  L = length(PC)
  # FDR level q
  q = 0.01
  w = which(sort(PC) < q * (1:L) / L)
  candidates.bh = order(PC)[w]
  filename <- paste0("Candidates_",name,"_q0.01.RDS")
  saveRDS(candidates.bh,filename)
  
  candidate_info<- data.frame(snp_number=candidates.bh,log10_p=-log10(PC)[candidates.bh],
                             snp_id=ann$snp_id[candidates.bh],effect=ann$effect[candidates.bh],
                             gene_ID=ann$ID[candidates.bh],gene_name=ann$gene_name[candidates.bh])
  write.csv(candidate_info,paste0("Candidates_info_",name,".csv"),row.names = F)
  
  # #Plotting
  # print("Plotting")
  # putative.snps <- which(PC_reduced<=tail(PC[candidates.bh],1))
  # 
  # filenamepdf <- paste0(name,"-LFMM_logp-value.pdf")
  # pdf(filenamepdf,width=10,height=5)
  # plot(-log10(PC_reduced$PC),
  #      pch = 19,
  #      cex = .2,
  #      xlab = "SNP", ylab = "-Log P",
  #      col = "grey")
  # points(putative.snps,
  #        -log10(PC_reduced$PC)[putative.snps],
  #        col = "blue",
  #        cex = .4)
  # dev.off()
}

######################### INVESTIGATE ANNOTATION ############################

candidates_PC1T <- readRDS("Candidates_PCT1_q0.01.RDS")
candidates_PC2T <- readRDS("Candidates_PCT2_q0.01.RDS")
candidates_PC3T <- readRDS("Candidates_PCT3_q0.01.RDS")
candidates_PC1P <- readRDS("Candidates_PCP1_q0.01.RDS")
candidates_PC2P <- readRDS("Candidates_PCP2_q0.01.RDS")
candidates_PC3P <- readRDS("Candidates_PCP3_q0.01.RDS")
env.variables <- list(candidates_PC1T,candidates_PC2T,candidates_PC3T,
                      candidates_PC1P,candidates_PC2P,candidates_PC3P)
names(env.variables) <- c("PC1T","PC2T","PC3T","PC1P","PC2P","PC3P")

# Loop through candidates for different variables

for(i in 1:length(env.variables)){
  candidates_PC <- env.variables[[i]]
  PC_snps <- ann[candidates_PC,]
  name <- names(env.variables)[i]
  print(name)
  
  # What effects do the candidates have?
  effect <- table(PC_snps$effect)
  effect <- effect[effect!=0]
  write.csv(effect,paste0(name,"_effect-table.csv"),row.names = F)
  pdf(paste0(name,"_effect-barplot.pdf"),width = 25,height = 10)
  barplot(effect,col=sample(rainbow(length(effect))))
  dev.off()
  
  dat <- as.data.frame(effect)
  pdf(paste0(name,"_effect-pizza.pdf"))
  p <- ggplot(dat, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() # remove background, grid, numeric labels
  print(p)
  dev.off()
  
  gene_id <- na.exclude(unique(PC_snps$ID))
  gene_name <- (unique(PC_snps$gene_name))
  gene_name <- gene_name[!gene_name==""][-1]
  snp_id <- PC_snps$snp_id
  write.csv(gene_id,paste0(name,"_gene_id.csv"),row.names = F)
  write.csv(gene_name,paste0(name,"_gene_name.csv"),row.names = F)
  write.csv(snp_id,paste0(name,"_snp_id.csv"),row.names = F)
}

# NOW FOR ALL CANDIDATES OF TEMP AND PRECIP

temp_candidates <- unique(union_all(candidates_PC1T,candidates_PC2T,candidates_PC3T))
prec_candidates <- unique(union_all(candidates_PC1P,candidates_PC2P,candidates_PC3P))
saveRDS(temp_candidates,"Candidates_TEMP_q0.01.RDS")
saveRDS(prec_candidates,"Candidates_PREC_q0.01.RDS")
length(intersect(temp_candidates,prec_candidates))
sets <- list(temp_candidates,prec_candidates)
names(sets) <- c("temperature","precipitation")

for(i in 1:length(sets)){
  candidates_PC <- sets[[i]]
  PC_snps <- ann[candidates_PC,]
  name <- names(sets)[i]
  print(name)
  
  # What effects do the candidates have?
  effect <- table(PC_snps$effect)
  effect <- effect[effect!=0]
  write.csv(effect,paste0(name,"_effect-table.csv"),row.names = F)
  pdf(paste0(name,"_effect-barplot.pdf"),width = 25,height = 10)
  barplot(effect,col=sample(rainbow(length(effect))))
  dev.off()
  
  dat <- as.data.frame(effect)
  pdf(paste0(name,"_effect-pizza.pdf"))
  p <- ggplot(dat, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() # remove background, grid, numeric labels
  print(p)
  dev.off()
  
  gene_id <- na.exclude(unique(PC_snps$ID))
  gene_name <- unique(PC_snps$gene_name)
  gene_name <- gene_name[!gene_name==""][-1]
  snp_id <- PC_snps$snp_id
  write.csv(gene_id,paste0(name,"_gene_id.csv"),row.names = F)
  write.csv(gene_name,paste0(name,"_gene_name.csv"),row.names = F)
  write.csv(snp_id,paste0(name,"_snp_id.csv"),row.names = F)
}

################ GO TERM ENRICHMENT ################

candidates_info_PC1T <- read.csv("Candidates_info_PCT1.csv")[,5]
candidates_info_PC2T <- read.csv("Candidates_info_PCT2.csv")[,5]
candidates_info_PC3T <- read.csv("Candidates_info_PCT3.csv")[,5]
candidates_info_PC1P <- read.csv("Candidates_info_PCP1.csv")[,5]
candidates_info_PC2P <- read.csv("Candidates_info_PCP2.csv")[,5]
candidates_info_PC3P <- read.csv("Candidates_info_PCP3.csv")[,5]
candidates_info_T <- list(candidates_info_PC1T,candidates_info_PC2T,candidates_info_PC3T)
candidates_info_P <- list(candidates_info_PC1P,candidates_info_PC2P,candidates_info_PC3P)

# # There is a lot going on here, but I basically had to write a clunky function to go over the list of 
# # gene IDs and get each individual ID (several of them were together because of intergenic regions)
# # (other had weird "path" in the middle of the name)
# get.gene.ids <- function(set=set){
#   complete.list.of.gene.ids <- numeric()
#   
#   for(i in 1:length(set)){
#     
#     candidates_info_PC <- set[[i]]
#     
#     all.genes <- numeric()
#     for(line in 1:length(candidates_info_PC)){
#       change.here <- unlist((strsplit(candidates_info_PC[line],"-")))
#       new.ids <- numeric()
#       for(item in change.here){
#         a <- unlist((strsplit(item,".path\\d")))
#         new.ids <- c(new.ids,a)
#       }
#       all.genes <- c(all.genes,new.ids)
#     }
#     # final lost of genes 
#     uniq.genes <- unique(all.genes)
#     complete.list.of.gene.ids <- c(complete.list.of.gene.ids,uniq.genes)
#   }
#   return(unique(complete.list.of.gene.ids))
# }
# 
# # This list includes all genes on either sides of an intergenic region
# temp.genes.ALL <- get.gene.ids(candidates_info_T)
# prec.genes.ALL <- get.gene.ids(candidates_info_P)
# write.csv(temp.genes.ALL,"TEMP.gene_IDS.csv",row.names = F)
# write.csv(prec.genes.ALL,"PREC.gene_IDS.csv",row.names = F)

temp.all.gene_names <- unique(union_all(read.csv("Candidates_info_PCT1.csv")[,6],read.csv("Candidates_info_PCT2.csv")[,6],
                                        read.csv("Candidates_info_PCT3.csv")[,6]))

prec.all.gene_names <- unique(union_all(read.csv("Candidates_info_PCP1.csv")[,6],read.csv("Candidates_info_PCP2.csv")[,6],
                                        read.csv("Candidates_info_PCP3.csv")[,6]))
write.csv(temp.all.gene_names,"temp.all.gene_names.csv",row.names = F)
write.csv(prec.all.gene_names,"prec.all.gene_names.csv",row.names = F)

go <- read.delim("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\PANNZER2\\PANNZER2 Gene Ontology\\GO.out.txt",sep="\t",header=TRUE)

####################################
#                                 
#
#            TEMPERATURE
#
#
###################################

# Here I'm considering only genes containing SNPs within or nearby
temp.gene <- unique(na.exclude(union_all(candidates_info_T[[1]],candidates_info_T[[2]],candidates_info_T[[3]])))
temp.genes <- temp.gene[nchar(temp.gene)==11] #this removed everything that does not have a single gene ID
write.csv(temp.genes,"temperature_gene_id.UNIQUE.csv",row.names = F)
temp.genes <- read.csv("temperature_gene_id.UNIQUE.csv")
temp.genes <- temp.genes$x

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
write.csv(BP,"GO.TEMP.BP.only_within.csv",row.names = F)

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
write.csv(fisher_test,"BP-temp-fisher_test.enrichment.csv",row.names = F)

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
prec.gene <- unique(na.exclude(union_all(candidates_info_P[[1]],candidates_info_P[[2]],candidates_info_P[[3]])))
prec.genes <- prec.gene[nchar(prec.gene)==11] #this removed everything that does not have a single gene ID
write.csv(prec.genes,"precipitation_gene_id.UNIQUE.csv",row.names = F)

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
write.csv(BP,"GO.PREC.BP.only_within.csv",row.names = F)

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


##########################################
#                                        #
#            Produce BED FILE            #
#                                        #
##########################################

temp_candidates <- readRDS("Candidates_TEMP_q0.01.RDS")
prec_candidates <- readRDS("Candidates_PREC_q0.01.RDS")

temp_len <- length(temp_candidates)
prec_len <- length(prec_candidates)

temp_bed <- data.frame(chrom=ann$chromosome.x[temp_candidates],chromStart=ann$position[temp_candidates],
                       chromEnd=ann$position[temp_candidates])
prec_bed <- data.frame(chrom=ann$chromosome.x[prec_candidates],chromStart=ann$position[prec_candidates],
                       chromEnd=ann$position[prec_candidates])
write.table(temp_bed,"temp_candidates.bed",sep = "\t",row.names = F)
write.table(prec_bed,"prec_candidates.bed",sep = "\t",row.names = F)

# Negative comparison (SNPs not associated to environment)

either_candidates <- unique(union(temp_candidates,prec_candidates))
total_number_of_SNPs <- nrow(ann)
seq <- 1:total_number_of_SNPs
NOTcandidate <- seq[-either_candidates]

set.seed(123458)
#temp_sample <- sample(NOTcandidate,temp_len)
temp_sample <- NOTcandidate
temp_NOTcandidates <- data.frame(chrom=ann$chromosome.x[temp_sample],chromStart=ann$position[temp_sample],
                                 chromEnd=ann$position[temp_sample])
#prec_sample <- sample(NOTcandidate,prec_len)
prec_sample <- NOTcandidate
prec_NOTcandidates <- data.frame(chrom=ann$chromosome.x[prec_sample],chromStart=ann$position[prec_sample],
                                 chromEnd=ann$position[prec_sample])
write.table(temp_NOTcandidates,"temp_NOTcandidates.bed",sep = "\t",row.names = F)
write.table(prec_NOTcandidates,"prec_NOTcandidates.bed",sep = "\t",row.names = F)
