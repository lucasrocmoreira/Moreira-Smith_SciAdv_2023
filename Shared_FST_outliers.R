downy <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Fst-outlier (ANGSD)\\Picoides pubescens\\Corrected SFS - v2\\Gene_id_all_outliers.csv",header=T)
hairy <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Fst-outlier (ANGSD)\\Picoides villosus\\SFS corrected - v2\\Gene_id_all_outliers.csv",header=T)

same_genes <- intersect(downy$x,hairy$x)

write.csv(same_genes,'Common_genes_fst_outlier_D&H.csv',row.names = F)

# Get Gene IDs

annotation <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\Annotation.csv",header=T)

annot.info <- data.frame()
for(i in same_genes){
  x <- annotation[grep(i,annotation$ID),]
  annot.info <- rbind(annot.info,x)
}
gene.id <- unique(annot.info$gene_name)[-1]
write.csv(gene.id,"Shared-gene_names.csv",row.names = F)

###############################################################

downy.unique.windows <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Fst-outlier (ANGSD)\\Picoides pubescens\\Corrected SFS - v2\\total.unique.windows.csv",header=T)
hairy.unique.windows <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Fst-outlier (ANGSD)\\Picoides villosus\\SFS corrected - v2\\total.unique.windows.csv",header=T)

downy.unique.windows$x <- paste0(downy.unique.windows$chr,"-",downy.unique.windows$midPos)
hairy.unique.windows$x <- paste0(hairy.unique.windows$chr,"-",hairy.unique.windows$midPos)

length(intersect(downy.unique.windows$x,hairy.unique.windows$x))
