##Collapser for DC annotation tables

##Purpose: if there is more than one segment per event, collapse into one event.
##If all the segments are ambiguous, set it as ambiguous overall.
##If at least one segment is incompatible, and none are compatible, set it as incompatible overall.
##If one segment is compatible, set it as compatible overall.


#set working directory to the source file location

library(plyr)

files = list.files(pattern="DC_annotated_table") #collecting the individual annotated tables for each genotype
files1 = length(files) # Count number of files
files2 = read.table(text = files, sep = "_", as.is = TRUE) #Split file names by "_" separator and create table "files2"

master <- NULL  #create a master table 
  
for (j in 1:files1) {
  data=NULL 
  data2=NULL
  data <- read.delim(files[j], sep = "\t", header=TRUE, fill=TRUE)
  filename <- files2[j,1]# save Meiosis identifier from files2 into a variable
  
  #test output
  #length(unique(data$EID))
  
  ##problem: event IDs are re-used across different meiosis##
  ##solution: extract the meiosis identifier and tag it to the front of the EID##
  for(q in 1:nrow(data)){
    tag=gsub('.*([0-9]+)','\\1',data[q,"Meiosis"])  # Extracts digit portion of string
    if(tag==0){tag<-10} #lazy way of fixing an issue with OMT10
    data[q,"EID"] <- paste(tag, data[q,"EID"], sep="")
  }

  for (e in unique(data$EID)){   
    subset <- data[which(data$EID == e),] #select each individual event
    
    if(nrow(subset)==1){  #if it only has 1 row (only 1 6:2 segment), just pass it unchanged
     data2 <- rbind(data2, subset)
    }
    
    if(nrow(subset)>1){ #if it has more than 1 6:2 segment....
      
        subset1 <- subset[grep("CC",subset$Classification),]  #make a subset for each classification
        subset2 <- subset[grep("CN",subset$Classification),]
        subset3 <- subset[grep("IN",subset$Classification),]
        subset4 <- subset[grep("IC",subset$Classification),]
        subset5 <- subset[grep("AN",subset$Classification),]
        subset6 <- subset[grep("AC",subset$Classification),]
        
        #CO:if there is at least 1 compatible 6:2, pass that line to the next table
        if(nrow(subset1)>0){                        
          data2 <- rbind(data2, subset1)
        }
        #NCO:if there is at least 1 compatible 6:2, pass that line to the next table
        if(nrow(subset2)>0){                        
          data2 <- rbind(data2, subset2)
        }
        #CO:if there aren't any compatible 6:2s but there are incompatible 6:2s, pass the first incompatible 6:2 to the next table
        if(nrow(subset3)>0 & nrow(subset1)==0 & nrow(subset2)==0){
          data2 <- rbind(data2, subset3[1,])
        }
        #NCO:if there aren't any compatible 6:2s but there are incompatible 6:2s, pass the first incompatible 6:2 to the next table
        if(nrow(subset4)>0 & nrow(subset1)==0 & nrow(subset2)==0){ 
          data2 <- rbind(data2, subset4[1,])
        }
        #CO:if there aren't any compatible 6:2s or incompatible 6:2s, but there are ambiguous 6:2s, pass the first ambiguous 6:2 to the next table
        if(nrow(subset5)>0 & nrow(subset1)==0 & nrow(subset2)==0 & nrow(subset3)==0 & nrow(subset4)==0){ 
          data2 <- rbind(data2, subset5[1,])
        }
        #NCO:if there aren't any compatible 6:2s or incompatible 6:2s, but there are ambiguous 6:2s, pass the first ambiguous 6:2 to the next table
        if(nrow(subset6)>0 & nrow(subset1)==0 & nrow(subset2)==0 & nrow(subset3)==0 & nrow(subset4)==0){
          data2 <- rbind(data2, subset6[1,])
        }
    }
  #write each individual table  
  write.table(data2,file=sprintf("%sCollapsed_DC_Annotated_Table.txt",filename),row.names=F, quote=F, sep='\t') 
  }
  master <- rbind.fill(master, data2) # Also gather each table into master table
}

write.table(master, file="Collapsed DC master table", col.names =TRUE, row.names = FALSE, quote = FALSE, sep="\t", append=TRUE)
