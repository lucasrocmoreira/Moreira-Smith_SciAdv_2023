## Written by Lucas R. Moreira
## Last updated 28 March 2020
## Usage: takes output from ANGSD fst windows analysis and produce outlier plots

#Before running, make sure to add a column name "fst" to the file:
#for i in *fst_slidingwindow; do sed -i 's/Nsites/Nsites\tfst/g' $i; done

#install.packages("qqman")
library(qqman)

############################################################
########          50kb window - 10kb steps          ########
############################################################

for(file in list.files(path=".",pattern=".fst_slidingwindow")){
  
  name <- unlist(strsplit(file, "\\."))[1]
  print(name)
  fst <- read.table(file, header=TRUE)
  fstsubset<-fst[complete.cases(fst),]  # make sure there is no missing data
  SNP<-c(1:(nrow(fstsubset)))           # give number to SNP-windows
  df<-cbind(SNP,fstsubset)            # create dataframe with SNP-window numbers

  #### Now we need to reorder and replace the chromosome names by numbers

  print("Reordering and replacing the chromosome names by numbers")
  chrOrder<-c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","LG2","LGE22","Z","w")
  df$chr<-factor(df$chr, levels=chrOrder)
  df <- df[order(df$chr),]
  
  num = 0
  elem = "nothing"
  chr = numeric()
  for(i in df$chr){
    if(i==elem){
      chr <- c(chr,num)
    }
    else{
      elem <- i
      num <-  num+1
      chr <- c(chr,num)
    }
  }
  
  df$CHROM <- chr                    # rename the scaffolds/chrom
  mydf <- df
  mydf <- mydf[mydf$chr!="w",]
  mydf <- mydf[mydf$chr!="Z",]
  
  # ##### Now we select only the windows with more than 20 SNPs
  # hist(mydf$Nsites,breaks = seq(min(mydf$Nsites),max(mydf$Nsites),
  #                          by=((max(mydf$Nsites) - min(mydf$Nsites))/(length(mydf$Nsites)-1))))
  # sd(mydf$Nsites)
  # mean(mydf$Nsites)+2*sd(mydf$Nsites)
  # mean(mydf$Nsites)-2*sd(mydf$Nsites)
  # 
  # new_mydf <- mydf[mydf$Nsites>33806,]
  
  ##### We then check the distribution of FSTs
  
  print("Checking distribution of FSTs")
  values <- mydf$fst
  mean_values <- mean(values)
  significance_sd <- mean_values+(5*sd(values))     # a suggested significance threshold (5 standard deviations above the mean)
  significance_top <- quantile(values, 0.99)                         # the top 1%
  pdf(paste0(name,".FST_histogram.pdf"),height=4,width=6)
  hist(values,breaks = seq(min(values),max(values),
                           by=((max(values) - min(values))/(length(values)-1))))
  abline(v=significance_sd,col="red")
  dev.off()
  
  ##### Now we make the Manhattan plot
  
  print("Making Manhattan plot")
  pdf(paste0(name,".50kb-10kbstep.pdf"),height=4,width=10)
  manhattan(mydf,chr="CHROM",bp="midPos",p="fst",snp="SNP",
            ylim = c(0, 0.8),col=c("blue4", "orange3"),cex = 0.6,logp=FALSE,suggestiveline = significance_sd,genomewideline = mean_values,
            ylab="Weir and Cockerham Fst",xlab="Chromosome",chrlabs =c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","LG2","LGE22"))
  dev.off()
  
  ##### Now we ask "What are the SNP-windows above the significance threshold?"
  
  print("Exporting outliers")
  outliers <- mydf[mydf$fst>significance_sd,]
  
  outliers_list <- fst[outliers$SNP,]
  write.csv(outliers_list,paste0(name,".outliers_5sd.csv"))

}

############# Zoom into chromosome 2 ################

fst <- read.table("PV-AK-SE.fst_slidingwindow", header=TRUE)

chr2 <- fst[fst$chr=="2",]

chr2_sub <- chr2[2000:8000,]
seq <- (2000:8000)*10

pdf("chr2-2k.8k.50kb-10kbstep.pdf",height=6,width=10)
plot(seq,chr2_sub$fst, xlab="Chromosome 2 position (kb)", ylab="FST",pch=16)
dev.off()

chr2_sub <- chr2[4500:6000,]
seq <- (4500:6000)*10
pdf("chr2-4.5k.6k.50kb-10kbstep.pdf",height=6,width=10)
plot(seq,chr2_sub$fst, xlab="Chromosome 2 position (kb)", ylab="FST",ylim=c(0,1),pch=16)
points(x=53221.613,y=1,col="red",pch=16)
points(x=53187.351,y=1,col="blue",pch=16)
dev.off()

######################################################

for(file in list.files(path=".",pattern=".fst_slidingwindow")){
  
  name <- unlist(strsplit(file, "\\."))[1]
  print(name)

}
