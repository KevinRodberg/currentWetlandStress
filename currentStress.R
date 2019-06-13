#==================================================================================================
# Y:/proj/CFWI_WetlandStress/Update2018/EvalCurrentStress/currentStress.R
#==================================================================================================
#==================================================================================================
# Probable Currently Stressed Class 3 Wetlands
#
# Created by Kevin A. Rodberg - June 2019
#
# Purpose: 
#
#==================================================================================================
#--
#   package management:
#     provide automated means for first time use of script to automatically 
#	  install any new packages required for this code, with library calls 
#	  wrapped in a for loop.
#--
pkgChecker <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}

list.of.packages <-c( "githubinstall",
                      "tcltk2","rModflow","future.apply","future","listenv","reshape")

suppressWarnings(pkgChecker(list.of.packages))

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if ("rModflow" %in% new.packages) devtools::install_github("KevinRodberg/rModflow")
lapply(list.of.packages,require, character.only=TRUE)

options(warn=-1)
#--------------------------------------------------------------------------------------------------
# Provides GUI to choose model  
# - may not be needed any long in this code since its used in P80headDifference.R
#--------------------------------------------------------------------------------------------------
source ("//ad.sfwmd.gov/dfsroot/data/wsd/SUP/devel/source/R/ReusableFunctions/tclFuncs.R")

MFmodel.Params <- defineMFmodel()
model <- chooseModel()
M <- as.data.frame(MFmodel.Params[model,])


readModflowTopo <- function(){
  topoDir <- "//whqhpc01p/hpcc_shared/dbandara/CFWI/ECFTX/Model/Transient/20190327_Calib_Final/model/dis_tr11"
  topoFile <-paste(topoDir,"TR9b__L1TOP.tr9b_nodry.dat",sep='/')
  topo<-scan(topoFile)
}
readHeadsbinByLay <- function(filPtr, selectLayer,maxSP) {
  bigVector <- vector('numeric')
  HeaderRead <- readHeadsHeader(filPtr)
  kntFloats <- HeaderRead$K * HeaderRead$NR * HeaderRead$NC
  Lay1floats <- HeaderRead$NR * HeaderRead$NC
  HeadBlock <- readBin(filPtr, double(), n = Lay1floats, size = 4)
  bigVector <- c(bigVector, HeadBlock[1:Lay1floats])
  i <- 1
  cat(paste("0%.."))
  SP_rng <- maxSP-1
  repeat {
    HeaderRead <- readHeadsHeader(filPtr)
    # Don't read past EOF
    if (length(HeaderRead) > 0) {
      if (HeaderRead$K == selectLayer) {
        i <- i + 1
        HeadBlock <-
          readBin(filPtr, double(), n = Lay1floats, size = 4)
        bigVector <- c(bigVector, HeadBlock[1:Lay1floats])
      } else {
        seek(filPtr, (Lay1floats * 4), origin = 'current')
      }
    }
    # don't read everything unless necessary
    if (length(HeaderRead) == 0) {
      cat('\n')
      break
    }
    
    if (HeaderRead$KPER > max(SP_rng)) {
      cat('\n')
      break
    }
    # Display % complete
    cat(paste('\r',format(as.numeric(HeaderRead$KPER) / max(SP_rng) * 100,digits = 2,nsmall = 2),"%"))
  }
  
  return(bigVector)
}
getBinaryFile <- function(M,byLayer){
  msg = paste('Identify Binary Heads file for :',M$MFmodel,'Layer',byLayer)
  winA <- tktoplevel()
  lbl.message <- tk2label(winA, text = msg, font = fontHeading)
  tkgrid(lbl.message, padx = 30)
  tkraise(winA)
  
  headsFile<-choose.files()
  
  tkdestroy(winA)
  
  if (length(headsFile) == 0) {
    exit("User cancelled HeadsFile choice")
  }
  return(headsFile)
}

processBinaryHeads<-function(headsFile,M,byLayer){

  to.read = file(headsFile, "rb")
  
  #===============================================
  # Estimate number of stress periods in Heads file
  #===============================================
  fileSz <- file.info(headsFile)$size
  TtlStrPd = fileSz / ( M$nlays * ((M$ncols * M$nrows * 4) + 44))
  
  #===============================================
  # Define range of Stress Periods to read
  #===============================================
  # SP_rng <- readRange()
  SP_rng <- seq(1:TtlStrPd)
  if (max(SP_rng) > TtlStrPd || min(SP_rng) < 1) {
    exit('Out of Range')
  }
  
  #===============================================
  # Retrieve Heads by Layer
  #===============================================
  to.read <- file(headsFile, "rb")
  selectLayer = byLayer
  maxSP <- as.integer(max(SP_rng))
  Layer <- readHeadsbinByLay(to.read, selectLayer, maxSP)
  close(to.read)

  HeadsMatrix<- array(Layer,c(M$ncols,M$nrows,maxSP))
  # PartialHeadsMatrix<- HeadsMatrix[,,122:133]

  plan(multiprocess)
  cat(paste('Calculating P80 heads from',maxSP,'stress periods with Layer',selectLayer,'Heads','\n'))
  xf <-future_apply (HeadsMatrix,MARGIN=c(1,2),FUN=stats::quantile,probs=c(.2),na.rm=T)
  # xf <-future_apply (PartialHeadsMatrix,MARGIN=c(1,2),FUN=stats::quantile,probs=c(.2),na.rm=T)
}
areaWeighted <- function(df,param){
  Summary <-merge(merge(aggregate(df$AcresMM, by=list(df$EMT_ID), FUN=sum),
                        aggregate(df$AcresMM*df[,c(param)], by=list(df$EMT_ID), FUN=sum),by="Group.1"),
                  aggregate(df[,c(param)], by=list(df$EMT_ID), FUN=mean), by="Group.1")
  Summary$newVal <- Summary[,3]/Summary[,2]
  names(Summary) <- c("EMT_ID","TotalArea",paste0(param,"XArea"),paste0("Avg",param),paste0("AreaWgt",param))
  return(Summary)
}
plotHistoDens <- function(fileName,C1,stress,phys){
  graphics.off()
  compare = melt(C1[C1$Phys==phys & C1$Status==stress,c('ObsTheta','SimTheta')])
  if(length(compare[compare$value==0,]$value)>0) {
    compare[compare$value==0,]$value<-NA
  }
  p <- ggplot(data=compare) +    
#    geom_histogram(aes(x=value,y=..density.., fill=variable),position="dodge",bins=15) +
    geom_histogram(aes(x=value, fill=variable),position="dodge",bins=15) +
    geom_density(aes(x=value,fill=variable),alpha=.2) +
    labs(title=paste(stress,phys),x = "Theta (Feet NAVD88)")
  ggsave(filename=fileName,width=10,height=6.66,units="in",dpi=300)
}

setwd('//ad.sfwmd.gov/dfsroot/data/wsd/SUP/proj/CFWI_WetlandStress/Update2018/EvalCurrentStress')
C1GrdPnts <- read.csv("Class1GridPnts.csv")
C1Pnts <- read.csv("Class1Pnts.csv")
C1ObsP80<- read.csv('AllP80.csv')
C2GrdPnts <- read.csv("Class2GridPnts.csv")

NotNeeded<-c("OBJECTID","ORIG_FID","Wetland__1", "Wetland_Ty","Column_")
C1GrdPnts[NotNeeded]<-NULL

names(C1GrdPnts) <-c("SEQNUM","EMT_ID", "Phys", "SHA", "ERE", "Urban_Density", 
                     "AcresMM", "XCOORD_UTM","YCOORD_UTM" )

NotNeeded <- c(  "OBJECTID","CFCA_ID_1","ACRES_2019","SiteName","Lake_Wetla","Lon","Lat" ,  
                 "FISSURE_V", "Subsiden_1","INVADING_V","EXP_ROOT_V","LEANFALL_V","DEAD_DYE_V",
                 "Success","Ridge","Acres","Kevin_s_Remark","Classification","Floodplain",
                 "Kym_s_Kris_s_Lisa_s_Comments","Row","Column_","ERE")

C2GrdPnts[NotNeeded]<-NULL

names(C2GrdPnts) <-c("EMT_ID","Phys","Stress", "Wetland_Type","Hydro_Class", 
                      "SEQNUM","AcresMM", "XCOORD_UTM","YCOORD_UTM" )

topo<-readModflowTopo()
heads1<-getBinaryFile(M,1)
p80HeadsLay1 %<-% processBinaryHeads(heads1,M,1)
heads3<-getBinaryFile(M,3)
p80HeadsLay3 %<-% processBinaryHeads(heads3,M,3)
# p80Heads<-xf
C1GrdPnts$topo <- topo[C1GrdPnts$SEQNUM]
C1GrdPnts$p80L1 <- p80HeadsLay1[C1GrdPnts$SEQNUM]
C1GrdPnts$p80L3 <- p80HeadsLay3[C1GrdPnts$SEQNUM]

C1Topo<-areaWeighted(C1GrdPnts,'topo')
C1P80L1<-areaWeighted(C1GrdPnts,'p80L1')
C1P80L3<-areaWeighted(C1GrdPnts,'p80L3')
C1<- merge(unique(C1GrdPnts[,c(2,3,4,5,6)]),
           merge(C1Topo[,c(1,5)],
                 merge(C1P80L1[,c(1,5)],C1P80L3[,c(1,5)])))

C2GrdPnts$topo <- topo[C2GrdPnts$SEQNUM]
C2GrdPnts$p80L1 <- p80HeadsLay1[C2GrdPnts$SEQNUM]
C2GrdPnts$p80L3 <- p80HeadsLay3[C2GrdPnts$SEQNUM]

C2Topo<-areaWeighted(C2GrdPnts,'topo')
C2P80L1<-areaWeighted(C2GrdPnts,'p80L1')
C2P80L3<-areaWeighted(C2GrdPnts,'p80L3')
C2<- merge(unique(C2GrdPnts[,c(1,2,3,4,5)]),
           merge(C2Topo[,c(1,5)],
                 merge(C2P80L1[,c(1,5)],C2P80L3[,c(1,5)])))

names(C2)[names(C2) == 'Stress'] <- 'Status'
C2$Status <-as.character(C2$Status)
C2$Status[which(C2$Status=="YES")] <- "Stressed"
C2$Status[which(C2$Status=="NO")] <- "Not Stressed"
C2$Phys <-as.character(C2$Phys)
C2$Phys[which(C2$Phys=="Plains")] <- "Plain"

names(C1Pnts)
C1<-merge(C1, C1Pnts[,c(3,7,14,15)], by.x="EMT_ID", by.y="CFCA_EMT_I")
names(C1)[names(C1) == 'Stress_Sta'] <- 'Status'
C1<-merge(C1, C1ObsP80[,c(2,6)])

names(C1ObsP80)
names(C1)
C1$ObsTheta<-C1$ERE-C1$`X2009.2017_P80`
C1$SimTheta<-NA
C1[C1$Phys=="Ridge",]$SimTheta<-
  C1[C1$Phys=="Ridge",]$AreaWgttopo-C1[C1$Phys=="Ridge",]$AreaWgtp80L3

C1[C1$Phys=="Plain",]$SimTheta<-
  C1[C1$Phys=="Plain",]$AreaWgttopo-C1[C1$Phys=="Plain",]$AreaWgtp80L1

C2$ObsTheta<-0
C2$SimTheta<-NA

C2[C2$Phys=="Ridge",]$SimTheta<-
  C2[C2$Phys=="Ridge",]$AreaWgttopo-C2[C2$Phys=="Ridge",]$AreaWgtp80L3

C2[C2$Phys=="Plain",]$SimTheta<-
  C2[C2$Phys=="Plain",]$AreaWgttopo-C2[C2$Phys=="Plain",]$AreaWgtp80L1

Cs<-rbind(C1[,names(C1) %in% names(C2)],
      C2[,names(C2) %in% names(C1)])

Cs$zScore = NA
for (phys in unique(C1$Phys)){
  for (stress in unique(C1$Status)){
    sdTheta <- 
      sd(Cs[Cs$Status==stress & Cs$Phys==phys,]$SimTheta, na.rm=TRUE)
    MeanTheta <- 
      mean(Cs[Cs$Status==stress & Cs$Phys==phys,]$SimTheta, na.rm=TRUE)

    Cs[Cs$Status==stress & Cs$Phys==phys,]$zScore <-
      (Cs[Cs$Status==stress & Cs$Phys==phys,]$SimTheta - MeanTheta)/sdTheta
    # Class 1 histograms
    # fileName = paste0('C:\\Users\\krodberg\\Desktop\\C1',stress,'_',phys,'_histo.png')
    # plotHistoDens(fileName,C1[C1$Status==stress &
    #                             C1$Phys==phys,],stress,phys)
    # # class 2 histograms
    # fileName = paste0('C:\\Users\\krodberg\\Desktop\\C2',stress,'_',phys,'_histo.png')
    # plotHistoDens(fileName,C2[C2$Status==stress &
    #                             C2$Phys==phys,],stress,phys)
    # # class 1 and class 2 combined histograms
    # fileName = paste0('C:\\Users\\krodberg\\Desktop\\Cs',stress,'_',phys,'_histo.png')
    # plotHistoDens(fileName,Cs[Cs$Status==stress &
    #                             Cs$Phys==phys,],stress,phys)
    # class 1 and class 2 combined histograms with outliers removed
    fileName = paste0('C:\\Users\\krodberg\\Desktop\\CsX',stress,'_',phys,'_histo.png')
    plotHistoDens(fileName,Cs[Cs$Status==stress & 
                                Cs$zScore < 2 & Cs$zScore > -2 &
                                Cs$Phys==phys,],stress,phys)
  }
}

