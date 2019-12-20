#Dedicated imager for computationally-identified double cut candidates

#Note1: this script is dependent on several files which would be produced by 
#the main recombination event mapping pipeline. See "Neale-Lab/OctadRecombinationMapping"
#Note2: in this script, it is necessary to edit the genotype name manually 
#in three places depending on which genotype you want to image


######################
master2<-read.delim("Final_Annotated_table_18-12-17.txt") 
#segs1<-read.delim("Candidate_listNFOMTkd.txt")
master2$Meiosis <- as.character(master2$Meiosis)
master2$debut <- as.numeric(master2$debut)
master2$fin <- as.numeric(master2$fin)
masterGR_noNA <- read.delim("MasterGRnoNATable_DC4")
masterGR <- read.delim("MasterGRTable_DC4")
masterGR$Var_len2 <- (masterGR$Var_len-1)
masterGR_noNA$Var_len2 <- (masterGR_noNA$Var_len-1)
#the variant lengths stored in this column are actually 1bp too long because they include the first base upstream. 
#We subtract 1 to make it more accurate.

master2 <- master2[complete.cases(master2$LCO),]

spo_file <- read.delim("WT_SCFullMap_Neeman.txt")
smooth_spo_file <- read.delim("Spo11Cer3PanHA100bp.txt")
TotalHits=sum(spo_file$Watson+spo_file$Crick)
spo_file$TotalHits=(spo_file$Watson+spo_file$Crick)
spo_file$TotalHpM=(spo_file$Watson+spo_file$Crick)/TotalHits*1000000 # Convert Pan data to Hits per million reads (HpM)
hotspots<- read.delim("WT_Hotspot.Table_Neeman.txt")
rmm_file<-read.delim("RMMSubSamp_simple")
AllElementsDUB = read.table("AllElementsDUB.txt", sep = "\t", header=TRUE) #Import datatable
rec_file<-read.delim("Rec8Marg.txt")
DC_midpoints <-read.delim("WT_DCFullMap_Neeman.txt")

orient<-read.delim("Orientation_table.txt", stringsAsFactors = FALSE)


##edit this as appropriate for your genotype###
segs<-read.delim("Candidate_list_OM.txt", stringsAsFactors = FALSE)  

# Subset the datatable
draw<-master2 # Temp holding table for events that will be drawn. 
  
#draw=subset(draw,threshold ==1500) # Threshold = 500,1500,5000 : Specifies which event calls to view. Event painting is controlled by the next variable: "Mainoverlaythreshold"
mainoverlaythreshold=1500 # The variable that specifies the threshold to use for main event painting/overlay

##edit this as appropriate for your genotype###
draw=subset(draw,Genotype =="OM" ) # Meiosis = OM2, OMT10, etc (refers to a specific meiosis)

###this uses the regular event table for drawing, but adds some extra information from the individual segment table####
##this is necessary as the segment table doesn't have different thresholds####

#remove extra variant lines
segs2 <- NULL
for(Q in 1:max(segs$SID)){ 
  subseg <- segs[which(segs$SID == Q),]
  if(nrow(subseg)>0){
     subseg2 <- subseg[1,]
     segs2 <- rbind(subseg2, segs2)
  }
}

###Make annotation table###
segs2$Classification <- 1
segs2$Notes <- ""
segs2 <- segs2[order(segs2$SID),]

##edit this as appropriate for your genotype###
write.table(segs2,file="OM_DC_annotation_table.txt", row.names=F, quote=F, sep='\t')

segs2$LocalEventHpM_NWT <- round(segs2$LocalEventHpM_NWT,2)
segs2$LocalEventHpM_DC <- round(segs2$LocalEventHpM_DC,2)
segs2$SegHpM_NWT <- round(segs2$SegHpM_NWT,2)
segs2$SegHpM_DC <- round(segs2$SegHpM_DC,2)
draw$segMinlength <- 0
draw$segMaxlength <-0
draw$EventHpMSC <-0
draw$EventHpMDC <-0
draw$SegHpMSC <- 0
draw$SegHpMDC <- 0
draw$HO_outcome <- 0
draw$SID <- 0
draw$len_cat <- 'N'
draw$seg_cat <- 'N'
draw$seg_start <- 0
draw$seg_stop <- 0

draw2 <- NULL
for(i in 1:nrow(segs2)){
  subset <- draw[which(draw$Meiosis==segs2[i,'Meiosis']),]
  subset2 <- subset[which(subset$id==segs2[i,'EID']),] ###this should introduce multiple lines of the same event if it has more than 1 candidate 6:2 segment
  subset2[1,'segMinlength'] <- segs2[i,"len_min"]
  subset2[1,'segMaxlength'] <- segs2[i,"len_max"]
  subset2[1,'EventHpMSC'] <- segs2[i,"LocalEventHpM_NWT"]
  subset2[1,'EventHpMDC'] <- segs2[i,"LocalEventHpM_DC"]
  subset2[1,'SegHpMSC'] <- segs2[i,"SegHpM_NWT"]
  subset2[1,'SegHpMDC'] <- segs2[i,"SegHpM_DC"]
  subset2[1,'HO_outcome'] <- segs2[i,"HO_outcome"]
  subset2[1,'SID'] <- segs2[i,"SID"]
  subset2[1,'len_cat'] <- segs2[i,"len_cat"]
  subset2[1,'seg_cat'] <- segs2[i,"seg_cat"]
  subset2[1,'seg_start'] <- segs2[i,"start5"]
  subset2[1,'seg_stop'] <- segs2[i,"stop3"]
  draw2 <- rbind(draw2,subset2)
}

draw2 <- draw2[order(draw2$SID),] 

#######################

ordre=1:8 # Order of plotting spores
extend=5000 # How many bp to extend the plotting either side of the event

# Automated PDF generation
wd = getwd(); out = paste(wd,"/","Output_Files","/","OMN_DC_Images_±", extend,"bp_",".pdf",sep=""); pdf(file=out, width=21,height=9);

#layout(matrix(c(0:16,17,17,18,18,0,0,19:34,35,35,36,36,0),44, 1, byrow = T))

# Replace funtion call with a loop that steps through the draw2 subtable

if(nrow(draw2)==0){print("no lines in draw2")}

if(nrow(draw2)>0){
  for (j in 1: nrow(draw2) ) {
    meiosis=draw2[j,"Meiosis"]
    chromo=draw2[j,"chr"]
    debut=draw2[j,"debut"]
    fin=draw2[j,"fin"]
    segdebut=draw2[j,"seg_start"]
    segfin=draw2[j,"seg_stop"]
    
    # How much to extend the plotting either side of the event
    #extend=5000 
    debut=debut-extend
    fin=fin+extend
    
    # graph_bin_annot4<-function(meiosis,chromo,debut,fin,ordre=1:8){
    
    # Sélection du format et des données
    # This sets up the graphing parameters: mfrow calculates number of rows of grpahics to be draw2n. mar sets the margins of the graph
    layout(matrix(c(0,1,1,1,2:17,18,18,19,19,20,20,21,22,22,23,23,24,24,24,0),35, 1, byrow = T))
    #layout.show(20)
    par(mar=c(0,6,0,2),oma = c(4,0,0,0),las=1) # Sets margins per graph and outside margins per grouped set
    
    # Subset the masterGR tables for the region of interest
    a<-subset(masterGR,Meiosis==meiosis & chr == chromo & (debut) <= pos_c & pos_c <= (fin)) 
    a1<-subset(masterGR_noNA,Meiosis==meiosis & chr == chromo & (debut) <= pos_c & pos_c <= (fin)) 
    a2<-subset(masterGR,Meiosis==meiosis & chr == chromo & (segdebut) <= pos_c & pos_c <= (segfin)) 
    
    a2<-subset(masterGR,Meiosis==meiosis & chr == chromo)  
    a2<-subset(a2,(segdebut) <= pos_c)
    a2<-subset(a2,(pos_c <= segfin))
    
    while (nrow(a1)<20){ # Check whether the region of interest contains sufficient SNPs to draw2 - this is important to prevent an error when a/a1 contain no rows due to rare zero SNP events (i.e. COs) that also lack SNPs in flanking thresholded region.
      debut=debut-1000
      fin=fin+1000
      a<-subset(masterGR,Meiosis==meiosis & chr == chromo & (debut) <= pos_c & pos_c <= (fin)) 
      a1<-subset(masterGR_noNA,Meiosis==meiosis & chr == chromo & (debut) <= pos_c & pos_c <= (fin)) 
    }
    
    if(nrow(a)>0){
      
      # Subset master based on Meiosis and Chromosome being draw2n
      e2=subset(master2, Meiosis==meiosis & chr == chromo)
      
      # Specify event colours based on type (CO or NCO) - fill in master table with this information
      for (i in 1:nrow(e2)){
        if(e2[i,"LCO"]==1){e2[i,"eventcolour"]="wheat"}
        if(e2[i,"LCO"]==0){e2[i,"eventcolour"]="thistle2"}
        if(e2[i,"LCO"]==2){e2[i,"eventcolour"]="wheat"} 
        if(e2[i,"LCO"]==3){e2[i,"eventcolour"]="wheat"}
        if(e2[i,"LCO"]=="U"){e2[i,"eventcolour"]="powderblue"}
      }
      
      #Subset e2 based on threshold values so that each can be plotted
      #e1000=subset(e2, Meiosis==meiosis & chr == chromo & threshold==1000)
     # e1500=subset(e2, Meiosis==meiosis & chr == chromo & threshold==1500)
      #e2000=subset(e2, Meiosis==meiosis & chr == chromo & threshold==2000)
      #mainoverlay=subset(e2,Meiosis==meiosis & chr == chromo)
      
      ######add SNP coordinates as text######
      plot(a$pos_c,a[,1],type="n",ylim=c(0,300),ylab=paste("SNPs"),cex.lab=1,font=1,xlab="",xlim=c(debut,fin), xaxt="n",yaxt="n", axes=FALSE, xaxs="i",yaxs="i")
      #type=n for no plotting. axes=FALSE gets rid of border. 
      text(a$pos_c, y=100, labels = a$pos_c, srt=90, pos=3, cex=1) 
      #cex=1 controls the size of the text; however, making it smaller doesn't really help with overlap
      arrows(a$pos_c,0,a$pos_c,50,col="black",length=0,lwd=1, code=3, lend=1)
      
      #######################################
      
      # Graphes des données brutes
      # Graphs of raw data
      for (i in (ordre+2)){  
        plot(a$pos_c,a[,i],type="l",lwd=2,col=2,cex.axis=0.8,ylim=c(0.6,300),yaxp=c(1,100,1),ylab=paste("S",i-2,sep=""),cex.lab=1,font=1,xaxs="i",yaxs="i",xlab="",xlim=c(debut,fin), xaxt="n", log="y")
        #ylab=paste("S",gsub("c","",names(fichier)[i]),"_c",chromo,sep="",collapse=NULL) # MCM's command
        
        # On impose xlim=c(debut,fin) pour caler avec les données de spo et rmm
        # We impose xlim=c(start,end) to wedge with spo data and rmm
        
        points(a$pos_c,a[,i],pch=4,col=2,cex=2)
        points(a$pos_c,a[,(i+8)],pch=4,col=4,cex=2)
        lines(a$pos_c,a[,(i+8)],col=4,lwd=2)
        points(a$pos_c,rep(100,length(a$pos_c)),pch=-1*a[,(i+16)]+17,cex=1.5,col=-2*a[,(i+16)]+4)
      } #End of i first loop
      
      # Graphes des données binarisées
      # Graphs of binary data
      position=c(85,95,85,95,85,95,85,95)
      haut<-0
      #############################
      orient1 <- orient[which(orient$Meiosis == meiosis),]  ###plotting strand orientation
      Values <- orient1[,(chromo+2)]
      #############################
      for (i in (ordre+2)){
        Number <- Values[(i-2)]
        plot(a$pos_c,rep(100,length(a$pos_c)),xaxt="n",yaxt="n",type="n",ylab=paste("S",i-2,sep=""),cex.lab=1,font=2,xaxs="i",yaxs="i",ylim=c(50,130),axes=FALSE,xlim=c(debut,fin))
        
        #############################
        par(new=TRUE)
        #plot(a[1,2],xaxt="n",yaxt="n",type="n",ylab=paste(Number),cex.lab=1,font=2,xaxs="i",yaxs="i",ylim=c(50,130),axes=FALSE,xlim=c(debut,fin))
        
        mtext(Number, side=2, line=1,las=2, col='black', cex.axis=0.7)
        
        #############################
        #TO DO~
        #values are sideways
        #replace NA with something else
        #Move values slightly to the right
        #Overplot the previous plot rather than creating a new line
        
        # On impose xlim=c(debut,fin) pour caler avec les données de spo et rmm
        alt<-ifelse(haut==0,120,120)
        bas<-ifelse(haut==0,60,60)
        arrows(a2$startSNP,90,a2$stopSNP,90,col="powderblue",length=0,lwd=30, lend=1) # Underlay a highlighting box across 6:2 segment region
        arrows(a$pos_c,position[i-2]-30,a$pos_c,position[i-2]+30,col=-2*a[,(i+16)]+4,length=0,lwd=1)
        
        #arrows(a$startSNP,position[i-2],a$stopSNP,position[i-2],col=-2*a[,(i+16)]+4,length=0,lwd=15, lend=1) # Add SNP start/stop boundary overlays
        arrows(a1$startSNP,position[i-2],a1$stopSNP,position[i-2],col=-2*a1[,(i+16)]+4,length=0,lwd=10, lend=1) # Add SNP start/stop boundary overlays
        
        haut<-1-haut
        
        
      } #End of i second loop
      
      # haut<-1-haut permet d'alterner les longueurs des traits entre haut = 0 et haut = 1
      # Top <-1 Top toggles the lengths of the lines between top = 0 and top = 1
      
      # Addition of plots that indicate the events specified by 1000 and 2000bp event-calling thresholds
      plot(a$pos_c,rep(100,length(a$pos_c)),xaxt="n",yaxt="n",type="n",ylab=paste("Events"),cex.lab=1.25,font=2,xaxs="i",yaxs="i",ylim=c(0,165),axes=FALSE,xlim=c(debut,fin))
      #arrows(e1000$debut,140,e1000$fin,140,col=e1000$eventcolour,length=0,lwd=12, code=3, lend=1); text(debut,140, "1.0 kb merge",pos=4) # Underlay a highlighting box across all event regions
      arrows(e2$debut,100,e2$fin,100,col=e2$eventcolour,length=0,lwd=12, code=3, lend=1); text(debut,100, "1.5 kb merge",pos=4) # Underlay a highlighting box across all event regions
      #arrows(e2000$debut,60,e2000$fin,60,col=e2000$eventcolour,length=0,lwd=12, code=3, lend=1); text(debut,60, "2.0 kb merge",pos=4) # Underlay a highlighting box across all event regions
      
      # Plot event midpoints
      #arrows(e1000$midpoint,130,e1000$midpoint,150,col="black",length=0,lwd=1, code=3, lend=1)
      arrows(e2$midpoint,90,e2$midpoint,110,col="black",length=0,lwd=1, code=3, lend=1)
      #arrows(e2000$midpoint,50,e2000$midpoint,70,col="black",length=0,lwd=1, code=3, lend=1)
      
      ####################################################################
      ###plot indels as triangles with length underneath###
      plot(a$pos_c,a[,1],type="n",ylim=c(0,300),ylab=paste("Indels"),cex.lab=1.25,font=2,xlab="",xlim=c(debut,fin), xaxt="n",yaxt="n", axes=FALSE, xaxs="i",yaxs="i")
      #type=n for no plotting. axes=FALSE gets rid of border. 
      aDEL=subset(a, type_k=="d")
      if(nrow(aDEL)>0){
        yDEL =rep(200, nrow(aDEL))
        points(aDEL$pos_c, yDEL, type = "p", pch=25, col="dimgrey", cex=2) 
        text(aDEL$pos_c, (yDEL-145), labels = aDEL$Var_len2, pos=3, cex=1) 
      }
      
      aINS=subset(a, type_k=="i")
      if(nrow(aINS)>0){
        yINS =rep(200, nrow(aINS))
        points(aINS$pos_c, yINS, type = "p", pch=24, col="dimgrey", cex=2) 
        text(aINS$pos_c, (yINS-145), labels = aINS$Var_len2, pos=3, cex=1) 
      }
      
      ####################################################################
      
      # Addition des données de Spo11 and RMM
      rbPal <- colorRampPalette(c('black','blue','red','pink')) # Creates a colour palette if chosing to plot Spo11 as a coloured histogram "h"
      spo_file$Col <- rbPal(10)[as.numeric(cut(spo_file$TotalHits,breaks = 10))]
      
      map<-subset(spo_file,Chr == chromo & debut <= Pos & Pos <= fin)
      map2<-subset(rmm_file,chr == chromo & debut <= pos & pos <= fin)
      #map3<-subset(red1_file,chr == chromo & debut <= pos & pos <= fin)
      #map4<-subset(nucleo_file,chr == chromo & debut <= pos & pos <= fin)
      map5<-subset(hotspots,Chr == chromo & debut <= Start & End <= fin) #Note: Pan file variables are called CHROM, HS_START, HS_END
      map6<-subset(rec_file,Chr == chromo & debut <= Position & Position <= fin)
     # map7<-subset(spo_file,Chr == chromo & debut <= Pos & Pos <= fin)
      #map7<-subset(smooth_spo_file,chr == chromo & debut <= pos & pos <= fin) # Smooth version not made for Neeman table yet
      map8 <- subset(DC_midpoints ,Chr == chromo & debut <= Mid & Mid <= fin)
      
      #plot(map$Pos,map$TotalHpM,type="l",xlim=c(debut,fin),ylim=c(10,max(map$TotalHpM)),yaxp=c(100,round(max(map$TotalHpM),-2),2),ylab="Spo11",bty="u",xaxs="i",lwd=2,col=map$Col,xaxt="n",cex.axis=0.9,cex.lab=1.25)
     # plot(map7$Pos,map7$TotalHpM,type="l",xlim=c(debut,fin),ylim=c(10,max(map7$TotalHpM)),yaxp=c(100,round(max(map7$TotalHpM),-2),2),ylab="Spo11",bty="u",xaxs="i",lwd=2,xaxt="n",cex.axis=0.9,cex.lab=1.25)
      
      ##############################
      # Hanninng window code:
      
      require("e1071") # This pacakge permits calculation of hanning function
      
      #Decompression code here. Because data is sparse, this is needed in order to populate all the missing rows prior to smoothing
      spo11.0=subset(spo_file,Chr == chromo & debut <= Pos & Pos <= fin) #Make a sub-table of the DSB data that only contains those rows where chr = 1 in range of interest
      spo11.1 <- data.frame(Chr=chromo, Pos=(debut:fin)) # Creates expanded dataframe with ONLY Chr and Pos locations
      spo11.1 <- merge(spo11.1,spo11.0, all=TRUE) # Merge expanded empty dataframe with compressed spo11.0 dataframe
      spo11.1[is.na(spo11.1)] <- 0 # Convert all NA values to zero
      
      win=101 #Smoothing Window length
      hw=hanning.window(win) #create hanning window (require package e1071 to be loaded)
      temp=c(rep(0,win),spo11.1$TotalHpM, rep(0,win)) # Extend by the length of the sliding window with zeros at both ends
      spo11.s=filter(temp,hw) # smooth the temp vector using the hann window and the filter function
      spo11.s=spo11.s[(win+1):(length(spo11.s)-win)] # trim smooth back to correct length
      
      spo11.1$Spo11smoothed=spo11.s #Write into starting spo11.0 tabel as new column
      
      
      #Decompression code here. Because data is sparse, this is needed in order to populate all the missing rows prior to smoothing
      #DC.0<- subset(DC_midpoints ,Chr == chromo & debut <= Mid & Mid <= fin)
      #DC.1 <- data.frame(Chr=chromo, Mid=(debut:fin)) # Creates expanded dataframe with ONLY Chr and Pos locations
      #DC.1 <- merge(DC.1,DC.0, all=TRUE) # Merge expanded empty dataframe with compressed DC.0 dataframe
      #DC.1[is.na(DC.1)] <- 0 # Convert all NA values to zero
      
      #win=101 #Smoothing Window length
      #hw=hanning.window(win) #create hanning window (require package e1071 to be loaded)
      #temp=c(rep(0,win),DC.1$TotalHpM, rep(0,win)) # Extend by the length of the sliding window with zeros at both ends
      #DC.s=filter(temp,hw) # smooth the temp vector using the hann window and the filter function
      #DC.s=DC.s[(win+1):(length(DC.s)-win)] # trim smooth back to correct length
      
      #DC.1$DCsmoothed=DC.s #Write into starting DC.0 tabel as new column
      
      plot(spo11.1$Pos,spo11.1$Spo11smoothed,type="l",xlim=c(debut,fin),ylim=c(0.1,max(spo11.1$Spo11smoothed)),yaxp=c(100,round(max(spo11.1$Spo11smoothed),-2),2),ylab=NA,bty="u",xaxs="i",lwd=2,xaxt="n",cex.axis=0.9,cex.lab=1.25,col='darkblue')
      mtext("Spo11", side=2, line=-4,las=2, col='darkblue', cex.axis=0.7)
      par(new=TRUE)
      #plot(DC.1$Mid,DC.1$DCsmoothed,type="l",xlim=c(debut,fin),ylim=c(0.1,max(DC.1$DCsmoothed)),yaxp=c(100,round(max(DC.1$DCsmoothed),-2),2),ylab=NA,yaxt="n",bty="u",xaxs="i",lwd=2,xaxt="n",cex.axis=0.9,cex.lab=1.25,col='darkred')
      plot(map8$Mid,map8$DCHpM,type="l",xlim=c(debut,fin),ylim=c(0.1,max(map8$DCHpM)),yaxt="n", ylab=NA, yaxp=c(100,round(max(map8$DCHpM),-2),2),bty="u",xaxs="i",lwd=2,xaxt="n",cex.axis=0.9,cex.lab=1.25,col='darkred')
      axis(4,cex.axis=0.9,cex.lab=1.25)
      mtext("DC", side=4, line=-2,las=2, col='darkred', cex.axis=0.7)
      
      
       ############################################
      plot(map5$MedianPoint,map5$Total_HpM,type="n",ylim=c(0,300),ylab="",cex.lab=1.25,font=2,xlab="",xlim=c(debut,fin), xaxt="n",yaxt="n", axes=FALSE, xaxs="i",yaxs="i")
      #type=n for no plotting. axes=FALSE gets rid of border.
      if(nrow(map5)>0){
        arrows(map5$Start,140,map5$End,140,col='lightsteelblue',length=0,lwd=12, code=3, lend=1); text="" # Underlay a highlighting box across all event regions
        text(map5$MedianPoint, y=0, labels = map5$TotalHpM, pos=3, cex=1) 
      }
      
    #  plot(map8$Mid,map8$DCHpM,type="n",ylim=c(0,300),ylab="",cex.lab=1.25,font=2,xlab="",xlim=c(debut,fin), xaxt="n",yaxt="n", axes=FALSE, xaxs="i",yaxs="i")
    #  #type=n for no plotting. axes=FALSE gets rid of border.
    #  if(nrow(map8)>0){
    #    arrows(map8$Mid,140,map8$Mid,140,col='lightsteelblue',length=0,lwd=12, code=3, lend=1); text="" # Underlay a highlighting box across all event regions
    #    text(map8$Mid, y=0, labels = map5$DCHpM, pos=3, cex=1) 
    #  }
      
      
      ##############################################################################################
      
      #plot(map4$pos,map4$nucleo,type="l",xlim=c(debut,fin),ylim=c(min(map4$nucleo),max(map4$nucleo)),yaxp=c(min(map4$nucleo),round(max(map4$nucleo),0),2),ylab="Nuc",bty="u",xaxs="i",lwd=4,col="grey",
      #     xaxt="n",cex.axis=0.9, cex.lab=1.5)
      #plot(map3$pos,map3$red1,type="l",xlim=c(debut,fin),ylim=c(0,max(map3$red1)),yaxp=c(0,round(max(map3$red1),0),2),ylab="Red1",bty="u",xaxs="i",lwd=4,col="red",
      #     xaxp=c((round(debut,-3)),(round(fin,-3)),(round(round(fin,-3)-round(debut,-3)) / 1000)),cex.axis=0.9, cex.lab=1.5)
      #plot(map3$pos,map3$red1,type="l",xlim=c(debut,fin),ylim=c(0,max(map3$red1)),yaxp=c(0,round(max(map3$red1),0),2),ylab="Red1",bty="u",xaxs="i",lwd=4,col="red",xaxt="n",cex.axis=0.9, cex.lab=1.5)
      
      ##########################################################################################################################################################
      #Now plot the gene datatrack
      #First subset the relevant data
      genes=AllElementsDUB #First make a copy of the ALLElements table
      genes=subset(genes,chr==chromo & start>(debut-10000) & stop<(fin+10000)) #Make a sub-table of ALLElements where chr = 1 and has limits just beyond plot range
      genes=subset(genes,type=="gene") #Make a sub-table of ALLElements
      #Now perform the plot
      plot(genes$start,genes$start, xaxt="n",yaxt="n",type="n", ylab=paste("Genes"),cex.lab=1.5,font=2, xlim=c(debut,fin), ylim=c(0,120),axes=FALSE) #set up empty plot
      # Following module draws arrows for each element
      xrange=fin-debut
      ahead=xrange/25 #make arrowhead length proportional to plot range
      ahead[(ahead>500)]=500 #limit max length to 500
      av=75 #arrow vertical location relative to plot dimensions
      ahw=15 #arrow/head width
      
      genesW=subset(genes,genename !="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
      if(nrow(genesW)>0){
        for (i in 1:nrow(genesW)){
          polygon(c(genesW[i,"start"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead,genesW[i,"stop"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead, genesW[i,"start"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="palegreen", border="palegreen4")
          text((genesW[i,"start"]+genesW[i,"stop"])/2,av, font=3, genesW[i,"genename"], cex=0.9) }
      }
      
      genesW=subset(genes,genename=="Dubious_ORF" & orientation =="+") #Make a sub-table of ALLElements
      if(nrow(genesW)>0){
        for (i in 1:nrow(genesW)){
          polygon(c(genesW[i,"start"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead,genesW[i,"stop"], genesW[i,"stop"]-ahead, genesW[i,"stop"]-ahead, genesW[i,"start"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="palegreen", border="palegreen4", lty=2)
          text((genesW[i,"start"]+genesW[i,"stop"])/2,av, font=3, genesW[i,"sysname"], cex=0.9) }
      }
      
      av=25 #arrow vertical location for Crick genes relative to plot dimensions
      genesC=subset(genes,genename !="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
      if(nrow(genesC)>0){
        for (i in 1:nrow(genesC)){
          polygon(c(genesC[i,"stop"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead,genesC[i,"start"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead, genesC[i,"stop"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="lightpink", border ="lightpink4")
          text((genesC[i,"start"]+genesC[i,"stop"])/2,av, font=3, genesC[i,"genename"], cex=0.9) }
      }
      
      genesC=subset(genes,genename=="Dubious_ORF" & orientation =="-") #Make a sub-table of ALLElements
      if(nrow(genesC)>0){
        for (i in 1:nrow(genesC)){
          polygon(c(genesC[i,"stop"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead,genesC[i,"start"], genesC[i,"start"]+ahead, genesC[i,"start"]+ahead, genesC[i,"stop"]),c(av+ahw,av+ahw,av+ahw+ahw,av,av-ahw-ahw,av-ahw, av-ahw), col="lightpink", border ="lightpink4", lty=2)
          text((genesC[i,"start"]+genesC[i,"stop"])/2,av, font=3, genesC[i,"sysname"], cex=0.9) }
      }
      ####Plot RMM track###
      plot(map2$pos,map2$rmm,type="l",xlim=c(debut,fin),ylim=c(1,max(map2$rmm)),yaxp=c(1,round(max(map2$rmm),0),2),ylab="RMM",bty="u",xaxs="i",lwd=4,col='darkgoldenrod',
           xaxp=c((round(debut,-3)),(round(fin,-3)),(round(round(fin,-3)-round(debut,-3)) / 500)),cex.axis=0.9, cex.lab=1.25)
      if(nrow(map6)>0){
        arrows(map6$Position,0,map6$Position,50,col="black",length=0,lwd=1, code=3, lend=1)
      }
      ###Make figure legend###
      
      title(xlab = paste(sep="","Red=S288c, Blue=SK1",
                         "\n",
                         meiosis, 
                         "  Chromosome=",chromo,
                         "  Event Min Len=", draw2[j,"len_min"], " bp",
                         "  Event Max Len=", draw2[j,"len_max"]," bp",
                         "  Seg Min Len=", draw2[j,"segMinlength"]," bp",
                         "  Seg Max Len=", draw2[j,"segMaxlength"]," bp",
                         " len_cat = ", draw2[j,"len_cat"],
                         " seg_cat = ", draw2[j,"seg_cat"],
                         " Event Id=", draw2[j,"id"],
                         "  Seg ID=", draw2[j,"SID"],
                         "\n",
                         "  Group=", draw2[j,"groupe"],
                         "  Class=", draw2[j,"classe"],
                         "  No. of COs=", draw2[j,"LCO"],
                         "  Variants=", draw2[j,"nb_snp"],
                         "  Event SC HpM=", draw2[j,"EventHpMSC"],
                         "  Event DC HpM=", draw2[j,"EventHpMDC"],
                         "  Segment SC HpM=", draw2[j,"SegHpMSC"],
                         "  Segment DC HpM=", draw2[j,"SegHpMDC"],
                         "  Hotspot Overlap=", draw2[j,"HO_outcome"],
                         "  Type=", draw2[j,"type"]),
            
            
            outer = T, line = 2, cex.lab=1.4) # Cex.lab controls the labelling size
    }
  } #End of j loop
}
dev.off()


