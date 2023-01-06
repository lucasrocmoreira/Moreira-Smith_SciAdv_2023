## Writen by Lucas R. Moreira
## Usage: Finds candidate SNP windows intersecting among methods

library(ggvenn)
library(data.table)
library(dplyr)

# Investigate how many windows are shared between ANGSDasso, PCAdapt, and H-scan

ANGSDasso <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\Window-based\\all.candidates.windows.csv")
PCAdapt <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\PCAdapt\\Manuscript_revision\\Picoides_pubescens\\all.candidates.windows.csv")
H_scan <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\H-scan\\Manuscript_revision\\Picoides_pubescens\\all_candidates_windows.csv")
H_scan.uniq <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\H-scan\\Manuscript_revision\\Picoides_pubescens\\all_unique_windows.csv")
total <- fread("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\Window-based\\PC1T.rol_win")

pdf("Venn_Diagram.pdf")
ggvenn(
  list("ANGSD-asso"=ANGSDasso$x,"PCAdapt"=PCAdapt$x,"H-scan"=H_scan$x), 
  fill_color = c("#CD534CFF", "#868686FF", "#EFC000FF"),
  text_size = 8,  stroke_size = 1, set_name_size = 8,
  show_percentage = F
)
dev.off()

# probability
three.way.overlap.count <- numeric()
for(i in 1:1000){
  print(i)
  window.set.1 <- sample(1:nrow(total),nrow(ANGSDasso))
  window.set.2 <- sample(1:nrow(total),nrow(PCAdapt))
  window.set.3 <- sample(1:nrow(total),nrow(H_scan))
  
  intersect <- intersect(intersect(window.set.1,window.set.2),window.set.3)
  three.way.overlap.count <- c(three.way.overlap.count,length(intersect))
}
hist(three.way.overlap.count,breaks=0:6,freq = T, right = F,xlab="Number of X windows in 3-way overlap",main=NULL)
probability <- sum(three.way.overlap.count>length(intersect(intersect(H_scan$x,PCAdapt$x),ANGSDasso$x)))+1/(1000+1)

# export three-way overlap
three.way.overlap <- sort(intersect(intersect(H_scan$x,PCAdapt$x),ANGSDasso$x))
write.csv(three.way.overlap,"three.way.overlap.csv",row.names = F)

# Annotation

intersection <- read.csv("three.way.overlap.csv")
ann <- fread("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\PP-SNP_annotation.csv",header=T,sep = ",")
sites <- fread("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\Window-based\\sites.tsv")
which.sites <- paste0(sites$Chromosome,"_",sites$Position)
ann.sites <- ann[ann$snp_id %in% which.sites,]

annotation.genes <- data.frame()
genes <- numeric()
for(i in 1:nrow(intersection)){
  print(paste(i,"of",nrow(intersection)))
  window.id <- intersection[i,]
  chrom <- unlist(strsplit(window.id,"_"))[1]
  start <- as.numeric(unlist(strsplit(window.id,"_"))[2])-24
  end <- as.numeric(unlist(strsplit(window.id,"_"))[2])+25
  ann.chr <- ann.sites[ann.sites$chromosome.x==chrom,]
  annotation <- ann.chr[start:end,]
  annotation.genes <- rbind(annotation.genes,annotation)
  annotation.genes <- annotation.genes %>% distinct()
  
  genes_ann <- annotation.genes$ID
  for(i in 1:length(genes_ann)){
    snp.gene <- genes_ann[i]
    if(grepl("-",snp.gene,fixed=T)){
      split.gene <- unlist(strsplit(snp.gene,"-"))
      genes <- unique(c(genes,split.gene))
    }else{
      genes <- unique(c(genes,snp.gene))
    }
  }
}  

write.csv(genes,"three.way_gene_id.csv",row.names = F)
write.csv(annotation.genes,"three.way_snp_annotation.csv",row.names = F)

gene.annotations <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\Annotation.csv")
genes <- read.csv("three.way_gene_id.csv")
colnames(genes) <- "ID"

x <- merge(genes,gene.annotations,by="ID")
write.csv(x,"three.way_gene_id.csv",row.names = F)
