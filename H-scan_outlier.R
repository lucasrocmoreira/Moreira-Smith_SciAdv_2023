library('dplyr')

hscanAK <- read.table("AK_all_chr_Hscan.out",header=T)
hscanE <- read.table("E_all_chr_Hscan.out",header=T)
hscanR <- read.table("R_all_chr_Hscan.out",header=T)
hscanNW <- read.table("NW_all_chr_Hscan.out",header=T)
hscan_list <- list(hscanAK,hscanE,hscanR,hscanNW)
names(hscan_list) <- c("AK","E","R","NW")

ann <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\PP-SNP_annotation.csv",header=T)

for(i in 1:length(hscan_list)){
  set <- hscan_list[[i]]
  name <- names(hscan_list)[i]
  print(name)
  
  significance_top <- quantile(set$H, 0.99) 
  
  print("Exporting outliers")
  outliers <- set[set$H>significance_top,]
  outliers$snp_id <- paste0(outliers$chr,"_",outliers$pos)
  
  candidate <- left_join(outliers, ann, by = c("snp_id" = "snp_id"))
  candidate <- candidate[order(candidate$H),]
  
  snp_candidate <- candidate$snp_id
  
  write.csv(candidate,paste0(name,".candidates_99quatile.csv"),row.names = F)
  write.csv(snp_candidate,paste0(name,".candidates_snps.csv"),row.names = F)
}

AK_candidates <- read.csv("AK.candidates_snps.csv")
E_candidates <- read.csv("E.candidates_snps.csv")
R_candidates <- read.csv("R.candidates_snps.csv")
NW_candidates <- read.csv("NW.candidates_snps.csv")

AK_unique <- setdiff(unlist(AK_candidates),unique(union_all(E_candidates$x,R_candidates$x,NW_candidates$x)))
E_unique <- setdiff(unlist(E_candidates),unique(union_all(AK_candidates$x,R_candidates$x,NW_candidates$x)))
R_unique <- setdiff(unlist(R_candidates),unique(union_all(E_candidates$x,AK_candidates$x,NW_candidates$x)))
NW_unique <- setdiff(unlist(NW_candidates),unique(union_all(E_candidates$x,R_candidates$x,AK_candidates$x)))

all_unique <- union_all(AK_unique,E_unique,R_unique,NW_unique)
write.csv(all_unique,"all_unique_snps.csv",row.names = F)

all_candidates <- unique(union_all(AK_candidates$x,E_candidates$x,R_candidates$x,NW_candidates$x))
write.csv(all_candidates,"all_candidates_snps.csv",row.names = F)
