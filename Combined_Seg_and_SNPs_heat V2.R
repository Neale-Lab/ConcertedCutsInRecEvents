##Use###
#Input: the master table made by the 'Combining_ALevents' script (non-split version)
#Output: 
#bad spores (non-octad) are removed, so no events from them are considered, but other spores from the same octad can be used
#Events are split into seperate segments and then into seperate SNPs
#Observed + expected Spo11 SC and DC hits are calculated for segments and SNPs

##Update Marg 02-07-18
#Removed bits dealing with individual SNPs (O/E and hotspot overlap)
#Added in more oligo tables
#added in hotspot overlap for max segment region

#Pan <- read.delim("FullMap.Cer3_Pan_HA_1_4h_c.txt")
NWT <- read.delim("WT_SCFullMap_Neeman.txt")
Ntel <- read.delim("tel1_SCFullMap_Neeman.txt")
Ntkd <-read.delim("tel1-kd_SCFullMap_Neeman.txt") 
SFA <-read.delim("F260A_SCFullMap_Sam.txt")  #from Sam Tischfield

DC_midpoints <-read.delim("WT_DCFullMap_Neeman.txt")

#setwd("/mnt/nfs2/gdsc/mc482/Event_Table_combining_and_imaging/ANEvents/Output_Files/")

#master<-read.delim("MasterEventTable_2017-07-24")
###Heat of all segments###
#master<-read.delim("Final_Annotated_table_18-12-17.txt")
#master<-read.delim("MasterNEventTable_woflank2017-09-15")
#master<-read.delim("MasterALEventTable_NoFlank2017-05-11")

master<-read.delim("masterAEventTable_16_8_18")
##Remove unwanted Genotypes###

#master <-master[which(master$Genotype != "OMT"),]
#master <-master[which(master$Genotype != "OMTkd"),]

#############Remove bad spores###############
meiosislist <- unique(master$Meiosis)
master2 <-NULL
for(i in 1:length(meiosislist)){  
  subset <- master[which(master$Meiosis==meiosislist[i]),]
  if(meiosislist[i]=='OM1'){subset <- subset[grep("7|8", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='OM5'){subset <- subset[grep("5|6", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='OM6'){subset <- subset[grep("1|2", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='OMT2'){subset <- subset[grep("3|4", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='OMT4'){subset <- subset[grep("1|2", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='OMT5'){subset <- subset[grep("5|6", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='OMT9'){subset <- subset[grep("3|4", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='OMTkd15'){subset <- subset[grep("5|6", subset$SpAff,invert=TRUE),]}
  ##Scott's samples
  if(meiosislist[i]=='WT1'){subset <- subset[grep("5|6", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='WT2'){subset <- subset[grep("1|2", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='WT5'){subset <- subset[grep("5|6", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='FA3'){subset <- subset[grep("7|8", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='FA5'){subset <- subset[grep("7|8", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='EA1'){subset <- subset[grep("3|4", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='EA2'){subset <- subset[grep("1|2", subset$SpAff,invert=TRUE),]}
  ##Bertrand's samples
  if(meiosislist[i]=='OMOM2'){subset <- subset[grep("1|2", subset$SpAff,invert=TRUE),]}
  if(meiosislist[i]=='OMHM2'){subset <- subset[grep("3|4", subset$SpAff,invert=TRUE),]}
  
  master2 <-rbind(master2, subset)
}


write.table(master2,file="List_of_all_clean_events.txt", row.names=F, quote=F, sep='\t')

#############Split events into segments###############
master <-master2
master <- master[grep("6|7|8", master$type), ] # Type contains 6,7, or 8
write.table(master,file="List_of_clean_DC_candidate_events_only.txt", row.names=F, quote=F, sep='\t')

#master <- master[which(master$nb_seg > 0),] #ones with '0' segs are actually complex events
master <- master[which(master$LCO >= 0),] #remove U/NA
master$seg_len_mid <- as.character(master$seg_len_mid)
master$type <-as.character(master$type)
master$seg_len_min <- as.character(master$seg_len_min)
master$seg_len_max <- as.character(master$seg_len_max)
master$seg_start5 <- as.character(master$seg_start5)
master$seg_start3 <- as.character(master$seg_start3)
master$seg_stop5 <- as.character(master$seg_stop5)
master$seg_stop3 <- as.character(master$seg_stop3)
master$seg_nb_SNP <- as.character(master$seg_nb_SNP)
master$SpAff <- as.character(master$SpAff)

#master <- master[which(master$id == 42),]
#master <- master[grep("4", master$type), ] # Type contains 6)


master1 <- master[which(master$Genotype == 'OM'),] 
#write.table(master1,file="OM_Chr_Removed.txt", row.names=F, quote=F, sep='\t')
master2 <- master[which(master$Genotype == 'OMT'),] 
master3 <- master[which(master$Genotype == 'TRM'),] 
master4 <- master[which(master$Genotype == 'TCMM'),]
master5 <- master[which(master$Genotype == 'OMTkd'),] 
master6 <- master[which(master$Genotype == 'OMN'),] 

master7 <- master[which(master$Genotype == 'WT'),] 
master8 <- master[which(master$Genotype == 'FA'),] 
master12 <- master[which(master$Genotype == 'EA'),] 

master9 <- master[which(master$Genotype == 'OMOM'),] 
master10 <- master[which(master$Genotype == 'OMHM'),] 

master11 <- master[which(master$Genotype == 'OEM'),] 

fileprocessor <-function(df, filename){
sub_master <- NULL
for (i in 1:nrow(df)){
  midlist <- df[i,'seg_len_mid']
  minlist <- df[i,'seg_len_min']
  maxlist <- df[i,'seg_len_max']
  typelist <- df[i,'type']
  start5list <-df[i,'seg_start5']
  start3list <-df[i,'seg_start3']
  stop5list <-df[i,'seg_stop5']
  stop3list <-df[i,'seg_stop3']
  Sp_Afflist <- df[i,'SpAff']
  seg_SNPslist <- df[i,'seg_nb_SNP']
  midlist <- strsplit(midlist, split= '_', fixed=TRUE)
  minlist <- strsplit(minlist, split= '_', fixed=TRUE)
  maxlist <- strsplit(maxlist, split= '_', fixed=TRUE)
  typelist <- strsplit(typelist, split= '_', fixed=TRUE)
  start5list <- strsplit(start5list, split= '_', fixed=TRUE)
  start3list <- strsplit(start3list, split= '_', fixed=TRUE)
  stop5list <- strsplit(stop5list, split= '_', fixed=TRUE)
  stop3list <- strsplit(stop3list, split= '_', fixed=TRUE)
  Sp_Afflist <- strsplit(Sp_Afflist, split= '_', fixed=TRUE)
  seg_SNPslist <- strsplit(seg_SNPslist, split= '_', fixed=TRUE)
  
  ##check if the typelist is too long - if so, remove last element
  #CO events cause a problem because they have more patterns than segments due to a 4:4 at the end
  df[i,"nb_seg"] <- lengths(maxlist) #complex events aren't given a nb_seg in MC's script
  if(lengths(typelist) > df[i,"nb_seg"]){
    x <- lengths(typelist)
    temp <- data.frame(typelist,typelist)
    remove <- nrow(temp)
    temp <- temp[-remove,]
    temp <- temp[,-1]
    typelist <- temp
    #a messy way of removing the last element from typelist when typelist is too long (as is the case for many types of CO)
  }
  
  df1 <- data.frame(midlist,minlist,maxlist,typelist,Sp_Afflist,seg_SNPslist)
  df2 <- data.frame(start5list, start3list, stop5list, stop3list)
  names(df1)[1] <- "len_mid"
  names(df1)[2] <- "len_min"
  names(df1)[3] <- "len_max"
  names(df1)[4] <- "type"
  names(df1)[5] <- "Sp_Aff"
  names(df1)[6] <- "seg_SNPs"
  names(df2)[1] <- "start5"
  names(df2)[2] <- "start3"
  names(df2)[3] <- "stop5"
  names(df2)[4] <- "stop3"
  df1$SID <- 1:nrow(df1)
  df2$SID <- 1:nrow(df2)
  if(nrow(df2)>1){ df2 <- head(df2,-1)}
  df1 <- merge(df1,df2,by="SID")
  df1$SID <- i
  df1$event_snp <- df[i,'nb_snp']
  df1$LocalEventHpM_PAN <- df[i,'LocalEventHpM_PAN']
  df1$LocalEventHpM_NWT <- df[i,'LocalEventHpM_NWT']
  df1$LocalEventHpM_Ntel <- df[i,'LocalEventHpM_Ntel']
  df1$LocalEventHpM_Ntkd <- df[i,'LocalEventHpM_Ntkd']
  df1$LocalEventHpM_SFA <- df[i,'LocalEventHpM_SFA']
  df1$LocalEventHpM_DC <- df[i,'LocalEventHpM_DC']
  #df$obs_exp <- df[i,24]
  df1$len_mid <- as.numeric(as.character(df1$len_mid))
  df1$total_nb_seg <- df[i,'nb_seg']
  df1$Event_len_max <- df[i,'len_max']
  df1$Event_len_mid <- df[i,'len_mid']
  df1$Event_len_min <- df[i,'len_min']
  temp <-df1[grep("6",df1$type),]
  df1$nb_62_seg <- nrow(temp)
  df1$LCO <- df[i,'LCO']
  df1$chr <-df[i,'chr']
  df1$Meiosis <- df[i,'Meiosis']
  df1$EID <- df[i,'id'] 
  df1$EventSpAff <- df[i,'SpAff'] 
  df1$LSB <- df[i,'LSB']
  df1$LNCO <- df[i,'LNCO']
  df1$groupe <-df[i,'groupe']
  df1$classe <-df[i,'classe']
  df1$Event_type <- df[i,'type']
  df1$debut <-  df[i,'debut']
  df1$fin <- df[i,'fin']
  
  sub_master <- rbind(sub_master,df1)

}
sub_master <- sub_master[order(sub_master$EID),]
sub_master$SID <- 1:nrow(sub_master)
write.table(sub_master,file=sprintf("%s_All_segments_list.txt",filename), row.names=F, quote=F, sep='\t')
}

fileprocessor(master1, "OM")
fileprocessor(master2, "OMT")
fileprocessor(master3, "TRM")
fileprocessor(master4, "TCMM")
fileprocessor(master5, "OMTkd")
fileprocessor(master6, "OMN")
fileprocessor(master7, "WT")
fileprocessor(master8, "FA")
fileprocessor(master9, "OMOM")
fileprocessor(master10, "OMHM")
fileprocessor(master11, "OEM")
fileprocessor(master12, "EA")

OM_master<-read.delim("OM_All_segments_list.txt")
OMT_master<-read.delim("OMT_All_segments_list.txt")
TRM_master<-read.delim("TRM_All_segments_list.txt")
TCMM_master<-read.delim("TCMM_All_segments_list.txt")
OMTkd_master<-read.delim("OMTkd_All_segments_list.txt")
OMN_master<-read.delim("OMN_All_segments_list.txt")

WT_master<-read.delim("WT_All_segments_list.txt")
FA_master<-read.delim("FA_All_segments_list.txt")
EA_master<-read.delim("EA_All_segments_list.txt")

OMOM_master<-read.delim("OMOM_All_segments_list.txt")
OMHM_master<-read.delim("OMHM_All_segments_list.txt")
OEM_master<-read.delim("OEM_All_segments_list.txt")

#####################################################################################################
#####################################################################################################
OM_master$Genotype <- 'OM'
OMT_master$Genotype <- 'OMT'
TRM_master$Genotype <- 'TRM'
TCMM_master$Genotype <- 'TCMM'
OMTkd_master$Genotype <- 'OMTkd'
OMN_master$Genotype <- 'OMN'

WT_master$Genotype <- 'WT'
FA_master$Genotype <- 'FA'
EA_master$Genotype <- 'EA'

OMOM_master$Genotype <- 'OMOM'
OMHM_master$Genotype <- 'OMHM'
OEM_master$Genotype <- 'OEM'

#####################################################################################################
#####################################################################################################

master_all <- rbind(OM_master, OMT_master, TRM_master, TCMM_master, OMTkd_master, OMN_master, WT_master, FA_master, EA_master, OMOM_master, OMHM_master, OEM_master)
master_all$start5 <- as.numeric(as.character(master_all$start5))
master_all$start3 <- as.numeric(as.character(master_all$start3))
master_all$stop5 <- as.numeric(as.character(master_all$stop5))
master_all$stop3 <- as.numeric(as.character(master_all$stop3))
master_all$len_max <- as.numeric(as.character(master_all$len_max))


##############Work out HpM for segments before retrieving SNPS so it can all be in one table##########
#############################################################
# Calculate Spo11 hits per segment
#Loop through events and calculate Pan HpM for each event, then add to table.
##Event HpM is calculated in the ALEvents Table Script###
##Now calculate Segment Pan HpM and compare it###

TotalHits=sum(NWT$Watson+NWT$Crick)
NWT$TotalHpM=(NWT$Watson+NWT$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)

TotalHits=sum(Ntel$Watson+Ntel$Crick)
Ntel$TotalHpM=(Ntel$Watson+Ntel$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)

TotalHits=sum(Ntkd$Watson+Ntkd$Crick)
Ntkd$TotalHpM=(Ntkd$Watson+Ntkd$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)

TotalHits=sum(SFA$Watson+SFA$Crick)
SFA$TotalHpM=(SFA$Watson+SFA$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)


master_all$SegHpM_NWT <- 0
master_all$SegHpM_DC <- 0
master_all$SegHpM_Ntel <- 0
master_all$SegHpM_Ntkd <- 0
master_all$SegHpM_SFA <- 0

masterSeg <-NULL
for(c in 1:16){
  NWTsubset <- NWT[which(NWT$Chr == c),]
  Ntelsubset <- Ntel[which(Ntel$Chr == c),]
  Ntkdsubset <- Ntkd[which(Ntkd$Chr == c),]
  SFAsubset <- SFA[which(SFA$Chr == c),]
  DCsubset <- DC_midpoints[which(DC_midpoints$Chr == c),]
  master_Segsubset <- master_all[which(master_all$chr == c),]
  
  for (i in 1:nrow(master_Segsubset)){
    NWT.1=subset(NWTsubset, Pos>=(master_Segsubset[i,"start5"]) & Pos <=(master_Segsubset[i,"stop3"]))
    master_Segsubset[i,"SegHpM_NWT"]=sum(NWT.1$TotalHpM)
    DC.1=subset(DCsubset, Mid>=(master_Segsubset[i,"start5"]) & Mid <=(master_Segsubset[i,"stop3"]))
    master_Segsubset[i,"SegHpM_DC"]=sum(DC.1$DCHpM)
    NTEL.1=subset(Ntelsubset, Pos>=(master_Segsubset[i,"start5"]) & Pos <=(master_Segsubset[i,"stop3"]))
    master_Segsubset[i,"SegHpM_Ntel"]=sum(NTEL.1$TotalHpM)
    NTKD.1=subset(Ntkdsubset, Pos>=(master_Segsubset[i,"start5"]) & Pos <=(master_Segsubset[i,"stop3"]))
    master_Segsubset[i,"SegHpM_Ntkd"]=sum(NTKD.1$TotalHpM)
    SFA.1=subset(SFAsubset, Pos>=(master_Segsubset[i,"start5"]) & Pos <=(master_Segsubset[i,"stop3"]))
    master_Segsubset[i,"SegHpM_SFA"]=sum(SFA.1$TotalHpM)

    
  }
  masterSeg <- rbind (masterSeg, master_Segsubset)
}


###How many hits would be expected to be in the segment region, based on the hits in the whole event###
masterSeg$exp_local_hits_seg_NWT <- 0
masterSeg$exp_local_hits_seg_Ntel <- 0
masterSeg$exp_local_hits_seg_Ntkd <- 0
masterSeg$exp_local_hits_seg_SFA <- 0
masterSeg$exp_local_hits_seg_DC <- 0

masterSeg$exp_local_hits_seg_NWT <- (masterSeg$LocalEventHpM_NWT/(masterSeg$Event_len_max))*masterSeg$len_max
masterSeg$exp_local_hits_seg_Ntel <- (masterSeg$LocalEventHpM_Ntel/(masterSeg$Event_len_max))*masterSeg$len_max
masterSeg$exp_local_hits_seg_Ntkd <- (masterSeg$LocalEventHpM_Ntkd/(masterSeg$Event_len_max))*masterSeg$len_max
masterSeg$exp_local_hits_seg_SFA <- (masterSeg$LocalEventHpM_SFA/(masterSeg$Event_len_max))*masterSeg$len_max
masterSeg$exp_local_hits_seg_DC <- (masterSeg$LocalEventHpM_DC/(masterSeg$Event_len_max))*masterSeg$len_max

#masterSeg$exp_local_hits_seg_NWT <- masterSeg$exp_local_hits_seg_NWT +0.000000000000000000001 #add a small constant to prevent infinite values
#masterSeg$exp_local_hits_seg_DC  <- masterSeg$exp_local_hits_seg_DC +0.000000000000000000001 #add a small constant to prevent infinite values


masterSeg$Seg_Obs_Exp_HpM_NWT=round(masterSeg$SegHpM_NWT/masterSeg$exp_local_hits_seg_NWT,4) # What is the fold difference between observed and expected association with Spo11 hits for each event?
masterSeg$Seg_Obs_Exp_HpM_Ntel=round(masterSeg$SegHpM_Ntel/masterSeg$exp_local_hits_seg_Ntel,4)
masterSeg$Seg_Obs_Exp_HpM_Ntkd=round(masterSeg$SegHpM_Ntkd/masterSeg$exp_local_hits_seg_Ntkd,4)
masterSeg$Seg_Obs_Exp_HpM_SFA=round(masterSeg$SegHpM_SFA/masterSeg$exp_local_hits_seg_SFA,4)
masterSeg$Seg_Obs_Exp_HpM_DC=round(masterSeg$SegHpM_DC/masterSeg$exp_local_hits_seg_DC,4) 

write.table(masterSeg,file="HpM_list.txt", row.names=F, quote=F, sep='\t')


######################################
##Determine overlap of max 6:2 segment region with a hotspot

hotspotsW <- read.delim("WT_Hotspot.Table.txt")  #WT from Neeman's paper
hotspotsW$Start <- as.numeric(hotspotsW$Start)
hotspotsW$End <- as.numeric(hotspotsW$End)
hotspotsF <- read.delim("hotspots_F260A.txt ") #from Sam Tischfield
hotspotsF$Start <- as.numeric(hotspotsF$Start)
hotspotsF$End <- as.numeric(hotspotsF$End)

loops<- read.delim("Rec8Marg.txt")

list<-read.delim("HpM_list.txt")

###this method of determining hotspot overlap asks if there is a hotspot within the maximum region of the 6:2 segment
#previously I was only asking if the actual called variants were in a hotspot.

list_hotspots <-NULL
meiosislist <- unique(list$Meiosis)

for(m in 1:length(meiosislist)){
  sublist <- list[which(list$Meiosis == meiosislist[m]),]
    if(nrow(sublist) > 0){
      sublist$boundary_cross<- 0
      sublist$hotspot_overlapW <- 0
      sublist$hotspot_boundary_crossW<- 0
      sublist$est_dist_DSBsW <- 0
      sublist$hotspot_overlapF <- 0
      sublist$hotspot_boundary_crossF<- 0
      sublist$est_dist_DSBsF <- 0}         #assign value of 0 to all segs  
    if(nrow(sublist) > 0){
      for(i in 1:nrow(sublist)){
        #does the segment cross a loop boundary (rec8 peak)?
        subloop <- loops[which(loops$Chr == sublist[i,"chr"] & loops$Position >= sublist[i,"start5"] & loops$Position <= sublist[i,"stop3"] ),] 
        if(nrow(subloop) >0){sublist[i,"boundary_cross"] <- 1}
        
        ###how many hotspots does each segment overlap? (wt)
        subhot1 <- hotspotsW[which(hotspotsW$Chr == sublist[i,"chr"] &  hotspotsW$Start <= sublist[i,"start5"] & hotspotsW$End >= sublist[i,"start5"] ),] #if hotspot overlaps start of event
        subhot2 <- hotspotsW[which(hotspotsW$Chr == sublist[i,"chr"] &  hotspotsW$Start <= sublist[i,"stop3"] & hotspotsW$End >= sublist[i,"stop3"] ),] #if hotspot overlaps end of event
        subhot3 <- hotspotsW[which(hotspotsW$Chr == sublist[i,"chr"] &  hotspotsW$Start >= sublist[i,"start5"] & hotspotsW$End <= sublist[i,"stop3"] ),] #if hotspot appears within event
        subhot4 <- hotspotsW[which(hotspotsW$Chr == sublist[i,"chr"] &  hotspotsW$Start <= sublist[i,"start5"] & hotspotsW$End >= sublist[i,"stop3"] ),] #if event appears within hotspot
        subhot <-rbind(subhot1, subhot2, subhot3, subhot4)
        ##remove duplicate entries as this can pick up the same hotspot multiple times
        
        subhot <- subset(subhot, !duplicated(subhot$TotalHpM)) 
        
        if(nrow(subhot) >0){sublist[i,"hotspot_overlapW"] <- nrow(subhot)}
        
        ##Are these hotspots ever in different loops?
        if(nrow(subhot)>1){  ##when there are >1 hotspots,
          ##are they ever seperated by a boundary? is there a rec8 peak between the end of one hotspot and the start of another? 
          subloop <- loops[which(loops$Chr == subhot[1,"Chr"] & loops$Position >= min(subhot$Start) & loops$Position <= max(subhot$End) ),] 
          
          sublist[i,"hotspot_boundary_crossW"] <- nrow(subloop)
          sublist[i,"est_dist_DSBsW"] <- max(subhot$Midpoint)-min(subhot$Midpoint)  #if we predict that the DSBs occur in different hotspots, what might be the distance between them?
        }
          ###how many hotspots does each segment overlap? (spo11-F260A)
          subhot1 <- hotspotsF[which(hotspotsF$chr == sublist[i,"chr"] &  hotspotsF$Start <= sublist[i,"start5"] & hotspotsF$End >= sublist[i,"start5"] ),] #if hotspot overlaps start of event
          subhot2 <- hotspotsF[which(hotspotsF$chr == sublist[i,"chr"] &  hotspotsF$Start <= sublist[i,"stop3"] & hotspotsF$End >= sublist[i,"stop3"] ),] #if hotspot overlaps end of event
          subhot3 <- hotspotsF[which(hotspotsF$chr == sublist[i,"chr"] &  hotspotsF$Start >= sublist[i,"start5"] & hotspotsF$End <= sublist[i,"stop3"] ),] #if hotspot appears within event
          subhot4 <- hotspotsF[which(hotspotsF$chr == sublist[i,"chr"] &  hotspotsF$Start <= sublist[i,"start5"] & hotspotsF$End >= sublist[i,"stop3"] ),] #if event appears within hotspot
          subhot <-rbind(subhot1, subhot2, subhot3, subhot4)
          ##remove duplicate entries as this can pick up the same hotspot multiple times
          
          subhot <- subset(subhot, !duplicated(subhot$avgWT_count)) 
          
          if(nrow(subhot) >0){sublist[i,"hotspot_overlapF"] <- nrow(subhot)}
          
          ##Are these hotspots ever in different loops?
          if(nrow(subhot)>1){  ##when there are >1 hotspots,
            ##are they ever seperated by a boundary? is there a rec8 peak between the end of one hotspot and the start of another? 
            subloop <- loops[which(loops$Chr == subhot[1,"chr"] & loops$Position >= min(subhot$Start) & loops$Position <= max(subhot$End) ),] 
            
            sublist[i,"hotspot_boundary_crossF"] <- nrow(subloop)
            sublist[i,"est_dist_DSBsF"] <- max(subhot$midpoint)-min(subhot$midpoint)  #if we predict that the DSBs occur in different hotspots, what might be the distance between them?
          }
      }
      list_hotspots <-rbind(list_hotspots, sublist)
  }
}

list_hotspots <- list_hotspots[order(list_hotspots$SID),]

write.table(list_hotspots,file="Hotspot_overlap.txt",row.names=F, quote=F, sep='\t')
        
######################
#marker for 6:2 segs or others
HOlist<-read.delim("Hotspot_overlap.txt")
HOlist$SegType <- 'Other'
HOlist$type <- as.character(HOlist$type)
for (i in 1:nrow(HOlist)){
if(HOlist[i,'type'] =='6:2'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='2:6'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='06:02'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='02:06'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='(6:2)'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='(2:6)'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='(2:6a)'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='(6:2a)'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='(2:6b)'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
if(HOlist[i,'type'] =='(6:2b)'){HOlist[i,'SegType'] <- '6:2/2:6 Segment'}
  
if(HOlist[i,'type'] =='(7:1)'){HOlist[i,'SegType'] <- '7:1/1:7 Segment'}
if(HOlist[i,'type'] =='(1:7)'){HOlist[i,'SegType'] <- '7:1/1:7 Segment'}
if(HOlist[i,'type'] =='(7:1a)'){HOlist[i,'SegType'] <- '7:1/1:7 Segment'}
if(HOlist[i,'type'] =='(1:7a)'){HOlist[i,'SegType'] <- '7:1/1:7 Segment'}
if(HOlist[i,'type'] =='(8:0)'){HOlist[i,'SegType'] <- '8:0/0:8 Segment'}
if(HOlist[i,'type'] =='(0:8)'){HOlist[i,'SegType'] <- '8:0/0:8 Segment'}
  
if(HOlist[i,'type'] =='(6:2i)'){HOlist[i,'SegType'] <- '6:2i/2:6i Segment'}
if(HOlist[i,'type'] =='(6:2ai)'){HOlist[i,'SegType'] <- '6:2i/2:6i Segment'}
if(HOlist[i,'type'] =='(6:2bi)'){HOlist[i,'SegType'] <- '6:2i/2:6i Segment'}
if(HOlist[i,'type'] =='(6:2ci)'){HOlist[i,'SegType'] <- '6:2i/2:6i Segment'}
if(HOlist[i,'type'] =='(2:6i)'){HOlist[i,'SegType'] <- '6:2i/2:6i Segment'}
if(HOlist[i,'type'] =='(2:6ai)'){HOlist[i,'SegType'] <- '6:2i/2:6i Segment'}
if(HOlist[i,'type'] =='(2:6bi)'){HOlist[i,'SegType'] <- '6:2i/2:6i Segment'}
if(HOlist[i,'type'] =='(2:6ci)'){HOlist[i,'SegType'] <- '6:2i/2:6i Segment'}

if(HOlist[i,'type'] =='4:4'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4)'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4a)'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4b)'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4CO)'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4aCO)'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4bCO)'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4cCO)'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4dCO)'){HOlist[i,'SegType'] <- '4:4 Segment'}
if(HOlist[i,'type'] =='(4:4eCO)'){HOlist[i,'SegType'] <- '4:4 Segment'}

if(HOlist[i,'type'] =='4:4ai'){HOlist[i,'SegType'] <- '4:4i Segment'}
if(HOlist[i,'type'] =='4:4bi'){HOlist[i,'SegType'] <- '4:4i Segment'}
if(HOlist[i,'type'] =='(4:4ai)'){HOlist[i,'SegType'] <- '4:4i Segment'}
if(HOlist[i,'type'] =='(4:4bi)'){HOlist[i,'SegType'] <- '4:4i Segment'}
if(HOlist[i,'type'] =='(4:4ci)'){HOlist[i,'SegType'] <- '4:4i Segment'}
if(HOlist[i,'type'] =='(4:4di)'){HOlist[i,'SegType'] <- '4:4i Segment'}
  

if(HOlist[i,'type'] =='5:3'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='05:03'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='5:3a'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='5:3b'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='5:3c'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='5:3d'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(5:3)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(5:3a)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(5:3b)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(5:3c)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(5:3d)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='3:5'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='03:05'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='3:5a'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='3:5b'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='3:5c'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='3:5d'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(3:5)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(3:5a)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(3:5b)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(3:5c)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
if(HOlist[i,'type'] =='(3:5d)'){HOlist[i,'SegType'] <- '5:3/3:5 Segment'}
  

}
            

write.table(HOlist,file="types_list.txt", row.names=F, quote=F, sep='\t')

#####Assign length category####
master<-read.delim("types_list.txt")

masterA <- master[which(master$SegType == '6:2/2:6 Segment'),]
masterB <- master[which(master$SegType == '7:1/1:7 Segment'),]
masterC <- master[which(master$SegType == '8:0/0:8 Segment'),]

master <- rbind(masterA, masterB, masterC)

master$seg_cat <- 'single'
for(c in 1:nrow(master)){
  if(master[c,"nb_62_seg"] >1){master[c, "seg_cat"] <- 'multi'}
}


#replace NA in obs/Exp with 0
is.na(master) <- do.call(cbind,lapply(master, is.infinite))
master[is.na(master)] <- 0


master$len_cat <- '>1000'

for (j in 1:nrow(master)){
  subset <- master[j,]
  if(subset$len_min <1001){ master[j, "len_cat"] <- '750_1000'}
  if(subset$len_min <751){ master[j, "len_cat"] <- '500_750'}
  if(subset$len_min <501){ master[j, "len_cat"] <- '300_500'}
  if(subset$len_min <301){ master[j, "len_cat"] <- '150_300'}
  if(subset$len_max >=30){ if(subset$len_min <151){ master[j, "len_cat"] <- '30_150'}}
  if(subset$len_max <30){ master[j, "len_cat"] <- '<30'}
}

master1 <- master[which(master$Genotype == 'OM'),] 
#write.table(master1,file="OM_Chr_Removed.txt", row.names=F, quote=F, sep='\t')
master2 <- master[which(master$Genotype == 'OMT'),] 
master3 <- master[which(master$Genotype == 'TRM'),] 
master4 <- master[which(master$Genotype == 'TCMM'),]
master5 <- master[which(master$Genotype == 'OMTkd'),] 
master6 <- master[which(master$Genotype == 'OMN'),] 

master7 <- master[which(master$Genotype == 'WT'),] 
master8 <- master[which(master$Genotype == 'FA'),] 
master12 <- master[which(master$Genotype == 'EA'),] 

master9 <- master[which(master$Genotype == 'OMOM'),] 
master10 <- master[which(master$Genotype == 'OMHM'),] 
master11 <- master[which(master$Genotype == 'OEM'),] 

write.table(master1,file="Candidate_list_OM.txt", row.names=F, quote=F, sep='\t')
write.table(master2,file="Candidate_list_OMT.txt", row.names=F, quote=F, sep='\t')
write.table(master3,file="Candidate_list_TRM.txt", row.names=F, quote=F, sep='\t')
write.table(master4,file="Candidate_list_TCMM.txt", row.names=F, quote=F, sep='\t')
write.table(master5,file="Candidate_list_OMTkd.txt", row.names=F, quote=F, sep='\t')
write.table(master6,file="Candidate_list_OMN.txt", row.names=F, quote=F, sep='\t')
write.table(master7,file="Candidate_list_WT.txt", row.names=F, quote=F, sep='\t')
write.table(master8,file="Candidate_list_FA.txt", row.names=F, quote=F, sep='\t')
write.table(master12,file="Candidate_list_EA.txt", row.names=F, quote=F, sep='\t')
write.table(master9,file="Candidate_list_OMOM.txt", row.names=F, quote=F, sep='\t')
write.table(master10,file="Candidate_list_OMHM.txt", row.names=F, quote=F, sep='\t')
write.table(master11,file="Candidate_list_OEM.txt", row.names=F, quote=F, sep='\t')

####make annotation tables###

genolist <- unique(master$Genotype)

for (B in 1:length(genolist)){  
  geno <- genolist[B]
  
  segs <- read.delim(file=sprintf("Candidate_list_%s.txt",geno),stringsAsFactors = FALSE)
  
  ###Make annotation table###
  segs$Classification <- 'AC'
  for(P in 1:nrow(segs)){
    if(segs[P,"LCO"] ==0){segs[P,"Classification"] <- 'AN'}
  }  
  
  #AC and AN stand for ambiguous CO and ambiguous NCO. You must manually annotate the events with the other categories
  
  segs$Notes <- ""
  segs <- segs[order(segs$SID),]
  
  #write.table(segs,file="WT_DC_annotation_table.txt", row.names=F, quote=F, sep='\t')
  write.table(segs,file=sprintf("%s_DC_annotation_table.txt",geno), row.names=F, quote=F, sep='\t')

}