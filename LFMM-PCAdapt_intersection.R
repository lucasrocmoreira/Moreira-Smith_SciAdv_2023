########################################
# Claculate number of orthologous SNPs #
########################################

PCAdapt <- readRDS("Outliers.RDS")

downy_temp <- readRDS("Candidates_TEMP_q0.01.RDS")
downy_prec <- readRDS("Candidates_PREC_q0.01.RDS")

intersect_temp <- intersect(PCAdapt,downy_temp)
intersect_prec <- intersect(PCAdapt,downy_prec)
saveRDS(intersect_temp,"intersect_TEMP.RDS")
saveRDS(intersect_prec,"intersect_PREC.RDS")

length(intersect_temp)
length(intersect_prec)

# NOW FOR ALL CANDIDATES OF TEMP AND PRECIP

intersect_temp <- readRDS("intersect_TEMP.RDS")
intersect_prec <- readRDS("intersect_PREC.RDS")
sets <- list(intersect_temp,intersect_prec)
names(sets) <- c("temperature","precipitation")

ann <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\PP-SNP_annotation.csv",header=T)
library(ggplot2)

for(i in 1:length(sets)){
  candidates_PC <- sets[[i]]
  PC_snps <- ann[candidates_PC,]
  name <- names(sets)[i]
  print(name)
  
  candidate_info<- data.frame(snp_number=candidates_PC,snp_id=ann$snp_id[candidates_PC],
                              effect=ann$effect[candidates_PC],gene_ID=ann$ID[candidates_PC],
                              gene_name=ann$gene_name[candidates_PC])
  write.csv(candidate_info,paste0("Candidates_info_",name,".csv"),row.names = F)
  
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

candidate_ID_temp <- read.csv("temperature_gene_id.csv")
candidate_ID_prec <- read.csv("precipitation_gene_id.csv")

go <- read.delim("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\PANNZER2\\PANNZER2 Gene Ontology\\GO.out.txt",sep="\t",header=TRUE)

####################################
#                                 
#
#            TEMPERATURE
#
#
###################################

temp.genes <- candidate_ID_temp$x[nchar(candidate_ID_temp$x)==11] #this removed everything that does not have a single gene ID
write.csv(temp.genes,"temperature_gene_id.UNIQUE.csv",row.names = F)

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
write.csv(fisher_test,"BP-fisher_test.enrichment.csv",row.names = F)

enrichment <- fisher_test[which(as.numeric(fisher_test$corrected.p)<0.05),]

write.csv(enrichment,"BP-enrichment.csv",row.names = F)

####################################
#                                 
#
#            PRECIPITATION
#
#
###################################

prec.genes <- candidate_ID_prec$x[nchar(candidate_ID_prec$x)==11] #this removed everything that does not have a single gene ID
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

#########################################
###     INTERSECTION WITH H-SCAN      ###
#########################################

hscan_snpid <- read.csv("H-scan_all_unique_snps.csv")[,1]

PCAdapt_snpid <- ann$snp_id[PCAdapt]

downy_temp_snpid <- ann$snp_id[downy_temp]
downy_prec_snpid <- ann$snp_id[downy_prec]

# LFMM vs H-scan

intersect_LFMM_H_temp <- intersect(downy_temp_snpid,hscan_snpid)
intersect_LFMM_H_prec <- intersect(downy_prec_snpid,hscan_snpid)

# PCAdapt vs H-scan

intersect_PCAdapt <- intersect(PCAdapt_snpid,hscan_snpid)

# Intersect all

intersect_all_temp <- intersect(intersect(PCAdapt_snpid,downy_temp_snpid),hscan_snpid)
intersect_all_prec <- intersect(intersect(PCAdapt_snpid,downy_prec_snpid),hscan_snpid)

########## What's going on in the 2-WAY SNPs? ##########

library('dplyr')

temp_ANN_PCAdapt_LFMM_SNPs <- left_join(data.frame(snp_id=ann$snp_id[intersect_temp]), ann, by = c("snp_id" = "snp_id"))
prec_ANN_PCAdapt_LFMM_SNPs <- left_join(data.frame(snp_id=ann$snp_id[intersect_prec]), ann, by = c("snp_id" = "snp_id"))
write.csv(temp_ANN_PCAdapt_LFMM_SNPs,"temp_ANN_PCAdapt_LFMM_SNPs.csv",row.names = F)
write.csv(prec_ANN_PCAdapt_LFMM_SNPs,"prec_ANN_PCAdapt_LFMM_SNPs.csv",row.names = F)


########## What's going on in the 3-WAY SNPs? ##########

library('dplyr')

temp_ANN_SNPs <- left_join(data.frame(snp_id=intersect_all_temp), ann, by = c("snp_id" = "snp_id"))
prec_ANN_SNPs <- left_join(data.frame(snp_id=intersect_all_prec), ann, by = c("snp_id" = "snp_id"))
write.csv(temp_ANN_SNPs,"temp_ANN_SNPs.csv",row.names = F)
write.csv(prec_ANN_SNPs,"prec_ANN_SNPs.csv",row.names = F)

