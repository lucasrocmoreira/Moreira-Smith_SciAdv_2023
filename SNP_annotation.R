#################################################
# Script to merge vcf data with gene annotation #
#################################################


readVariableWidthFile <- function(filePath){
  con <-file(filePath)
  lines<- readLines(con)
  close(con)
  slines <- strsplit(lines," ")
  colCount <- max(unlist(lapply(slines, length)))
  
  FileContent <- read.table(filePath,
                            header = FALSE,
                            fill = TRUE)
  return(FileContent)
}

# Downy

snp <- readVariableWidthFile("PP-SNP-only_simplified.ann.snpsummary")

library('dplyr')
snp <- snp %>% mutate_all(na_if,"")
colnames(snp) <- c("snp_number","chromosome","position","snp_id","effect","impact","ID")
snp$ID <- as.character(snp$ID)

# read annotation

annotation <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\Annotation.csv",header=T)
annotation$ID <- as.character(annotation$ID)
annotation_clean <- annotation %>% distinct(ID, .keep_all = TRUE)
# merge tables

#snp_ann <- merge(snp,annotation,by.y="ID")

snp_ann <- left_join(snp, annotation_clean, by = c("ID" = "ID"))
snp_ann <- snp_ann[,-c(8,9)]

write.csv(snp_ann,"PP-SNP_annotation.csv",row.names = F)

# Hairy

snp <- readVariableWidthFile("PV-SNP-only_simplified.ann.snpsummary")

library('dplyr')
snp <- snp %>% mutate_all(na_if,"")
colnames(snp) <- c("snp_number","chromosome","position","snp_id","effect","impact","ID")
snp$ID <- as.character(snp$ID)

# read annotation

annotation <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\Annotation.csv",header=T)
annotation$ID <- as.character(annotation$ID)
annotation_clean <- annotation %>% distinct(ID, .keep_all = TRUE)
# merge tables

#snp_ann <- merge(snp,annotation,by.y="ID")

snp_ann <- left_join(snp, annotation_clean, by = c("ID" = "ID"))
snp_ann <- snp_ann[,-c(8,9)]

write.csv(snp_ann,"PV-SNP_annotation.csv",row.names = F)
