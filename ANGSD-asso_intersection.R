#######################################################################
#                                                                     #
#    Script to find overlapping candidate genes between two species   #
#                   Written by Lucas R. Moreira                       #
#                                                                     #
#######################################################################

#################### GENES ####################

PCT_downy_candidate <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\Window-based\\temp.all.gene_ids.csv")
PCP_downy_candidate <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\Window-based\\prec.all.gene_ids.csv")

PCT_hairy_candidate <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_villosus\\Window-based\\temp.all.gene_ids.csv")
PCP_hairy_candidate <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_villosus\\Window-based\\prec.all.gene_ids.csv")

temperature_genes <- na.exclude(intersect(PCT_downy_candidate$x,PCT_hairy_candidate$x))
length(temperature_genes)

precipitation_genes <- na.exclude(intersect(PCP_downy_candidate$x,PCP_hairy_candidate$x))
length(precipitation_genes)

write.csv(temperature_genes,"parallel-temperature_genes.csv",row.names = F)
write.csv(precipitation_genes,"parallel-precipitation_genes.csv",row.names = F)

# PARALLELISM #
# GENES #

### Temperature ###

pdf("Gene-temp-hypergenometric.pdf",width = 5,height = 5)
plot(dhyper(0:10,nrow(PCT_downy_candidate),
            15137-nrow(PCT_downy_candidate),nrow(PCT_hairy_candidate)),type="l",xlab = "Number of Overlapping Genes", 
     ylab ="Probability of X Overlapping Genes",lwd=3)
abline(v=length(temperature_genes)-1,col="red",lwd=3, lty=2)
p <- phyper(length(temperature_genes)-1,nrow(PCT_downy_candidate),
            15137-nrow(PCT_downy_candidate),nrow(PCT_hairy_candidate),lower.tail = F)
text(x=10,y=0.45,labels=paste0("p = ",format(round(p, 9), nsmall = 3)))
dev.off()

# Alternatively
overlap.count.temp <- numeric()
for(i in 1:1000){
  print(i)
  window.set.1 <- sample(1:15137,nrow(PCT_downy_candidate))
  window.set.2 <- sample(1:15137,nrow(PCT_hairy_candidate))
  
  intersect <- intersect(window.set.1,window.set.2)
  overlap.count.temp <- c(overlap.count.temp,length(intersect))
}
hist(overlap.count.temp,breaks=0:12,freq = T, right = F,xlab="Number of X windows in 3-way overlap",main=NULL)
probability <- sum(overlap.count.temp>length(temperature_genes))+1/(1000+1)

### Precipitation ###

pdf("Gene-prec-hypergenometric.pdf",width = 5,height = 5)
plot(dhyper(0:8,nrow(PCP_downy_candidate),
            15137-nrow(PCP_downy_candidate),nrow(PCP_hairy_candidate)),type="l",xlab = "Number of Overlapping Genes", 
     ylab ="Probability of X Overlapping Genes",lwd=3)
abline(v=length(precipitation_genes)-1,col="red",lwd=3, lty=2)
p <- phyper(length(temperature_genes)-1,nrow(PCP_downy_candidate),
            15137-nrow(PCP_downy_candidate),nrow(PCP_hairy_candidate),lower.tail = F)
text(x=8.5,y=0.65,labels=paste0("p = ",format(round(p, 10), nsmall = 3)))
dev.off()

# Alternatively
overlap.count.prec <- numeric()
for(i in 1:1000){
  print(i)
  window.set.1 <- sample(1:15137,nrow(PCP_downy_candidate))
  window.set.2 <- sample(1:15137,nrow(PCP_hairy_candidate))
  
  intersect <- intersect(window.set.1,window.set.2)
  overlap.count.prec <- c(overlap.count.prec,length(intersect))
}
hist(overlap.count.prec,breaks=0:12,freq = T, right = F,xlab="Number of X windows in 3-way overlap",main=NULL)
probability <- sum(overlap.count.prec>length(precipitation_genes))+1/(1000+1)

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
