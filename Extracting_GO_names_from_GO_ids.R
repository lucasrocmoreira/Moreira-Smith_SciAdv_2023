# # load the GO library
# library(GO.db)
# 
# # extract a named vector of all terms
# goterms <- Term(GOTERM)
# 
# # work with it in R, or export it to a file
# write.table(goterms, sep="\t", file="goterms.txt")

GOTERMS <- read.table("goterms.txt",row.names = NULL)

# Extracting GO term from description

getGO_ID <- function(file=file){
  desciptions <- file$Functional.Category
  
  GO_info <- data.frame()
  
  for(i in 1:length(desciptions)){
    print(i)
    description2 <- substr(tolower(desciptions[i]),1,nchar(desciptions[i])-1)
    term <- GOTERMS[grep(paste0("^",description2,"$"),GOTERMS$x,ignore.case=TRUE),][,1] # "^" and "$" guarantees an exact match
    if(length(term)==0){
      term <- "NA"
    }
    y <- cbind(GO_id = term,Description=description2)
    GO_info <- rbind(GO_info,y)
    print(term)
  }
  
  return(GO_info)
}

Downy_finch <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Fst-outlier (ANGSD)\\Picoides pubescens\\Corrected SFS - v2\\ShinyGO-enrichment_ZEBRAFINCH_005.csv")
Downy_human <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Fst-outlier (ANGSD)\\Picoides pubescens\\Corrected SFS - v2\\ShinyGO-enrichment_HUMAN_005.csv")

GO_results1 <- getGO_ID(file=Downy_finch)
write.csv(GO_results1,"ShinyGO-enrichment_ZEBRAFINCH_005-GOIDS.csv",row.names = F)

GO_results2 <- getGO_ID(file=Downy_human)
write.csv(GO_results2,"ShinyGO-enrichment_HUMAN_005-GOIDS.csv",row.names = F)

Hairy_finch <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Fst-outlier (ANGSD)\\Picoides villosus\\SFS corrected - v2\\ShinyGO-enrichment_ZEBRAFINCH_005.csv")
Hairy_human <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Fst-outlier (ANGSD)\\Picoides villosus\\SFS corrected - v2\\ShinyGO-enrichment_HUMAN_005.csv")

GO_results3 <- getGO_ID(file=Hairy_finch)
write.csv(GO_results3,"ShinyGO-enrichment_ZEBRAFINCH_005-GOIDS.csv",row.names = F)

GO_results4 <- getGO_ID(file=Hairy_human)
write.csv(GO_results4,"ShinyGO-enrichment_HUMAN_005-GOIDS.csv",row.names = F)

files <- list.files("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\ShinyGO",pattern="ShinyGo_ZEBRAFINCH_FDR0_05")
outpath <- "C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\LFMM2\\ShinyGO"

for(i in files){
  complete.path <- paste0(outpath,"\\",i)
  read.this <- read.csv(complete.path)
  x <- getGO_ID(file=read.this)
  write.csv(x,paste0(outpath,"\\GO_",i),row.names = F)
}
