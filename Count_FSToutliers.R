#install.packages('reshape')
library(reshape)

total.count <- numeric()
position <- data.frame()
for(file in list.files(path="./FST comparisons/",pattern="..outliers_5sd")){
  print(file)
  name <- unlist(strsplit(file, "[.]"))[1]
  pop1 <- unlist(strsplit(name, "-"))[2]
  pop2 <- unlist(strsplit(name, "-"))[3]
  
  table <- read.csv(paste0("./FST comparisons/",file))
  pos <- table[,3:4]
  position <- rbind(position,pos)
  count = 1
  for(line in 2:nrow(table)){
    print(line)
    if(table$midPos[line]!=(table$midPos[line-1]+10000)){
      count = count + 1
    }
  }

  absolute.count = length(table$midPos)

  row.final <- cbind(pop1=pop1,pop2=pop2,countiguous=count,absolute=absolute.count)
  total.count <- rbind(total.count,row.final)
  
}
write.csv(total.count,"total.count.csv",row.names = F)

total.unique.windows <- position[!duplicated(position),]
write.csv(total.unique.windows,"total.unique.windows.csv",row.names = F)
