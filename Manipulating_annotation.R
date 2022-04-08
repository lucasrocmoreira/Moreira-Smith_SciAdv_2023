annotation <- read.table("Picoides_pubescens.gff")
genes <- annotation[annotation$V3=="mRNA",]
ID <- as.character(genes$V9)


library(stringr)
split_1 <- str_split(ID, pattern = fixed("="), simplify = TRUE, n = 4)
gene_number <- str_replace_all(split_1[,2], pattern = ";Source", replacement = "")
ENSEMBL <- str_replace_all(split_1[,3], pattern = ";Function", replacement = "")
gene_name <- str_replace_all(split_1[,4], pattern = ";", replacement = "")
gene_name <- str_replace_all(gene_name, pattern = "\"", replacement = "")

annotation_final <- data.frame(ID=gene_number,ENSEMBL=ENSEMBL,gene_name=gene_name)

###################################################################################

# Let's get the location of the genes

location <- read.table("Picoides_pubescens.pseudogenome.cds.gff.summary")
colnames(location) <- c("chromosome","first","last","ID")
split_2 <- str_split(location$ID , pattern = fixed("."), simplify = TRUE, n = 2)[,1]
split_3 <- str_split(split_2 , pattern = fixed("="), simplify = TRUE, n = 2)[,2]
location$ID <- split_3

###################################################################################

# Merge both tables

annotation_complete <- merge(location,annotation_final,by="ID")
annotation_complete <- annotation_complete[order(annotation_complete$chromosome,annotation_complete$first),]
write.csv(annotation_complete,"Annotation.csv")
