
## Script to test whether parallel genes are enriched for large coding sequences
## Written by Lucas R. Moreira

# All genes

all.genes <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\Annotation\\Annotation.csv")

# ANGSD-asso candidate genes

PP.temp <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\Window-based\\temp.all.gene_ids.csv")
PP.prec <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_pubescens\\Window-based\\prec.all.gene_ids.csv")
PV.temp <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_villosus\\Window-based\\temp.all.gene_ids.csv")
PV.prec <- read.csv("C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Chapter 3\\ANGSDasso\\Picoides_villosus\\Window-based\\prec.all.gene_ids.csv")

all.candidates <- unique(c(PP.temp$x,PP.prec$x,PV.temp$x,PV.prec$x))

# Classify top 25% genes as large

hist(all.genes$length)
quantile(all.genes$length,c(0.25,0.5,0.75,0.9))

classifier <- function(x){
  if(x>21408){
    return("L")
  } else{
    return("N")
  }
}

classification <- unlist(lapply(all.genes$length,classifier))

all.genes$classification <- classification

# count how many candidate genes are large

candidate.class <- all.genes[all.genes$ID %in% all.candidates,]$classification
sum(candidate.class=="L")
sum(candidate.class!="L")

# count how many total genes are large

sum(all.genes$classification=="L")
sum(all.genes$classification!="L")

sum(all.genes$classification=="L")-sum(candidate.class=="L")
sum(all.genes$classification!="L")-sum(candidate.class!="L")

# Fisher's Exact Test

dat <- data.frame(
  "large" = c(121,3953),
  "non-large" = c(403,11832),
  row.names = c("candidate","non-candidate"),
  stringsAsFactors = F
)
colnames(dat) <- c("Large", "Not large")
dat

mosaicplot(dat,
           main = "Mosaic plot",
           color = TRUE
)

fisher.test(dat,alternative = "greater")

#### Classifying top 50% as large

classifier <- function(x){
  if(x>6716){
    return("L")
  } else{
    return("N")
  }
}

classification <- unlist(lapply(all.genes$length,classifier))

all.genes$classification <- classification

# count how many candidate genes are large

candidate.class <- all.genes[all.genes$ID %in% all.candidates,]$classification
sum(candidate.class=="L")
sum(candidate.class!="L")

# count how many total genes are large

sum(all.genes$classification=="L")
sum(all.genes$classification!="L")

sum(all.genes$classification=="L")-sum(candidate.class=="L")
sum(all.genes$classification!="L")-sum(candidate.class!="L")

# Fisher's Exact Test

dat <- data.frame(
  "large" = c(202,7952),
  "non-large" = c(322,7833),
  row.names = c("candidate","non-candidate"),
  stringsAsFactors = F
)
colnames(dat) <- c("Large", "Not large")
dat

mosaicplot(dat,
           main = "Mosaic plot",
           color = TRUE
)

fisher.test(dat,alternative = "greater")
