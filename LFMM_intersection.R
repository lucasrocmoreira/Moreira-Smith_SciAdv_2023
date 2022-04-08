########################################
# Calculate number of orthologous SNPs #
########################################

hairy <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\PV-SNP_annotation.csv",header=T)

downy <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\PP-SNP_annotation.csv",header=T)

intersect <- intersect(hairy$snp_id,downy$snp_id)

length(intersect)

# add q0_001\\ to make it more stringent
PCT_downy <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\Downy\\Complete SNP set\\temperature_snp_id.csv",header=T)
length(intersect(intersect,PCT_downy$x))
PCP_downy <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\Downy\\Complete SNP set\\precipitation_snp_id.csv",header=T)
length(intersect(intersect,PCP_downy$x))

PCT_hairy <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\Hairy\\Complete SNP set\\temperature_snp_id.csv",header=T)
length(intersect(intersect,PCT_hairy$x))
PCP_hairy <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\Hairy\\Complete SNP set\\precipitation_snp_id.csv",header=T)
length(intersect(intersect,PCP_hairy$x))

temperature <- intersect(PCT_downy$x,PCT_hairy$x)
precipitation <- intersect(PCP_downy$x,PCP_hairy$x)

#################### GENES ####################


PCT_downy_candidate <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\Downy\\Complete SNP set\\temperature_gene_id.UNIQUE.csv")
PCP_downy_candidate <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\Downy\\Complete SNP set\\precipitation_gene_id.UNIQUE.csv")

PCT_hairy_candidate <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\Hairy\\Complete SNP set\\temperature_gene_id.UNIQUE.csv")
PCP_hairy_candidate <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\Hairy\\Complete SNP set\\precipitation_gene_id.UNIQUE.csv")

temperature_genes <- na.exclude(intersect(PCT_downy_candidate$x,PCT_hairy_candidate$x))
length(temperature_genes)

precipitation_genes <- na.exclude(intersect(PCP_downy_candidate$x,PCP_hairy_candidate$x))
length(precipitation_genes)

write.csv(temperature_genes,"parallel-temperature_genes.csv",row.names = F)
write.csv(precipitation_genes,"parallel-precipitation_genes.csv",row.names = F)

# PARALLELISM #
# SNP #

# 1-phyper(overlap(center),left_side,total,right_side)

### Temperature ###

pdf("SNP-temp-hypergenometric.pdf",width = 5,height = 5)
plot(dhyper(0:5,48,149423,110),type="l",xlab = "Number of Overlapping Outliers", 
     ylab ="Probability of X Overlapping Outliers",lwd=3)
abline(v=0,col="red",lwd=3, lty=2)
dev.off()
1-phyper(0,48,149423,110)

### Precipitation ###

pdf("SNP-prec-hypergenometric.pdf",width = 5,height = 5)
plot(dhyper(0:10,1322,149423,242),type="l",xlab = "Number of Overlapping Outliers", 
     ylab ="Probability of X Overlapping Outliers",lwd=3)
abline(v=4,col="red",lwd=3, lty=2)
dev.off()
1-phyper(4,1322,149423,242)

# GENES #

### Temperature ###

pdf("Gene-temp-hypergenometric.pdf",width = 5,height = 5)
plot(dhyper(0:250,714,15407,1204),type="l",xlab = "Number of Overlapping Genes", 
     ylab ="Probability of X Overlapping Genes",lwd=3)
abline(v=216,col="red",lwd=3, lty=2)
dev.off()
1-phyper(216,714,15407,1204)

### Precipitation ###

pdf("Gene-prec-hypergenometric.pdf",width = 5,height = 5)
plot(dhyper(0:2000,7218,15407,404),type="l",xlab = "Number of Overlapping Genes", 
     ylab ="Probability of X Overlapping Genes",lwd=3)
abline(v=1957,col="red",lwd=3, lty=2)
dev.off()
1-phyper(1957,7218,15407,404)

###################################
#                                 #
#         GO ENRICHMENT           #
#                                 #
###################################

###################################
#           TEMPERATURE           #
###################################

go <- read.delim("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\PANNZER2\\PANNZER2 Gene Ontology\\GO.out.txt",sep="\t",header=TRUE)
temp.genes <- temperature_genes
  
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
write.csv(fisher_test,"BP-temp-fisher_test.enrichment.csv",row.names = F)

enrichment <- fisher_test[which(as.numeric(fisher_test$corrected.p)<0.05),]

write.csv(enrichment,"BP-temp-enrichment.csv",row.names = F)

###################################
#          PRECIPITATION          #
###################################

prec.genes <- precipitation_genes

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

######################################### GET GENE NAMES ##############################################

library('dplyr')

ann <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\Annotation.csv",header=T)

temp <- read.csv("parallel-temperature_genes.csv")
prec <- read.csv("parallel-precipitation_genes.csv")

temp_ann <- left_join(temp, ann, by = c("x" = "ID"))
write.csv(temp_ann,"parallel-temperature_genesID.csv",row.names = F)
prec_ann <- left_join(prec, ann, by = c("x" = "ID"))
write.csv(prec_ann,"parallel-precipitation_genesID.csv",row.names = F)
