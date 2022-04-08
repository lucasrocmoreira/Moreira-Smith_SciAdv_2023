library(dplyr)

annotation <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\Annotation.csv",header=T)

go <- read.delim("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\PANNZER2\\PANNZER2 Gene Ontology\\GO.out.txt",sep="\t",header=TRUE)

############################ Get gene list from annotation ############################

genes <- data.frame()

for(file in list.files(path=".",pattern=".outliers_5sd.csv")){
  
  #file="PP-AK-MW.outliers_5sd.csv"
  print(file)
  a <- read.csv(file,header=T)
  a <- a[a$chr!="Z",]
  
  candidates <- data.frame()
  for(line in 1:nrow(a)){
    print(paste0(line," of ",nrow(a)))
    single_line <- unlist(strsplit(as.character(a$region[line]),split="[(:punct:),]"))
    first <- as.integer(single_line[5])
    last <- as.integer(single_line[6])
    chromosome <- a$chr[line]
    
    sub_annot <- annotation[annotation$chromosome==as.character(chromosome),]
    which_genes <- sub_annot[sub_annot$first>=first&sub_annot$first<=last,]
    
    candidates <- rbind(candidates,which_genes)
    
    genes <- rbind(genes,which_genes)
    genes <- genes %>% distinct()
  }
  final.candidates <- candidates[!duplicated(candidates),]
  file.name <- paste0(unlist(strsplit(as.character(file),split="[.]"))[1],"-candidates.csv")
  write.csv(final.candidates,file.name,row.names = F)
  
  # Now let's find the gene ontologies for the candidate genes
  gene.ontologies <- data.frame()
  for(i in 1:length(final.candidates$ID)){
    candidate <- final.candidates$ID[i]
    print(candidate)
    go_set <- go[grep(candidate,go$qpid),]
    ##go_set
    #go_CC <- go_set[grep("CC",go_set$ontology),]
    #CC <- rbind(CC,go_CC)
    #go_MF <- go_set[grep("MF",go_set$ontology),]
    #MF <- rbind(MF,go_MF)
    go_BP <- go_set[grep("BP",go_set$ontology),]
    gene.ontologies <- rbind(gene.ontologies,go_BP)
  }
  
  file.name2 <- paste0(unlist(strsplit(as.character(file),split="[.]"))[1],"-candidatesGO.csv")
  write.csv(gene.ontologies,file.name2,row.names = F)
}

uniq_genes <- genes[order(genes$chromosome,genes$first),]

write.csv(uniq_genes,"Annotation_of_all_outliers.csv",row.names = F)

gene_list <- unique(na.exclude(as.character(uniq_genes$ID)))
write.csv(gene_list,"Gene_id_all_outliers.csv",row.names = F)

###################### GO TERM ENRICHMENT ######################

library(GOfuncR)

Outlier_gene_ids = gene_list
input_hyper_outlier = data.frame(Outlier_gene_ids, is_candidate=1)
#res_hyper_outlier = go_enrich(input_hyper_outlier, orgDb='org.Gg.eg.db') #chicken
res_hyper_outlier = go_enrich(input_hyper_outlier) # humans does better
## first element of go_enrich result has the stats
stats = res_hyper_outlier[[1]]
## top-GO categories
head(stats)
## top GO-categories per domain
by(stats, stats$ontology, head, n=30)
## all valid input genes
head(res_hyper_outlier[[2]])
nrow(res_hyper_outlier[[2]])

##################### Get GO ID from terms extracted from ShinyGo (http://bioinformatics.sdstate.edu/go/) #####################

shinygo <- read.csv("ShinyGo_outlier_enrichment_BiologicalProcess_above005.csv",header=T)

library(GO.db)

goterms <- Term(GOTERM)
all_terms <- as.data.frame(goterms)

go_id <- numeric()
for(i in shinygo$Functional.Category){
  id <- rownames(all_terms)[which(all_terms==i)]
  if(length(id)<1){id <- NA}
  go_id <- c(go_id,id)
}

shinygo$GO_ID <- go_id

write.csv(shinygo,"ShinyGo_outlier_enrichment_BiologicalProcess_above005_withGOID.csv",row.names = F)
