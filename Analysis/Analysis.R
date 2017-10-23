
#Data
load("HTAData_RASA_PSR_LiverVSMuscle_AllPSR.RData")
load("HTA_ExonLevelSummarized_rma_RASA.RData")
load("HTA_GeneLevelSummarized_rma_RASA.RData")
colnames(HTA_GeneLevelSummarized_rma_RASA)=c("GeneID","Liver_A","Liver_B","Liver_C","Muscle_A","Muscle_B","Muscle_C")
load("REIDS_RankedProbesets.RData")
load("HTAData_RASA.rda")
for(i in 3:ncol(HTAData_RASA)){
	HTAData_RASA[,i]=as.numeric(as.character(HTAData_RASA[,i]))
}

PlotFunction<-function(GeneID=NULL,ExonID=NULL,Data,GDS_Output,GeneLevelData=NULL,ExonLevelData=NULL,FIRMA_Data=NULL,groups,ylabel=NULL,plottype="new",location){
	message("The gene and exon level data will be log2 transformed")
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			dev.off()
		}
	}
	
	
	if(is.null(GeneID) | is.null(ExonID)){
		stop("no GeneID and/or ExonID specified")
	}
	
	# Array Scores
	ArrayScores=GDS_Output[[which(names(GDS_Output)==GeneID)]]$arrayScore # the arrayscores for the entire gene
	
	# Exon Scores
#	ExonScores=GDS_Output[[which(names(GDS_Output)==GeneID)]]$exonScore
	
	# Exon Level Data
#	ExonLevel=ExonLevelData[which(ExonLevelData$GeneID==GeneID),]
#	ExonLevel=ExonLevel[,-c(1,2)]
#	rownames(ExonLevel)=ExonLevel[,2]
#	
#	ExonLevel=ExonLevel[which(rownames(ExonLevel%in%rownames(ArrayScores)),]
#	message("The Exon Level Data is log2 transformed")
#	ExonLevel=log2(ExonLevel)
#	ExonLevel=ExonLevel[,unlist(Groups)]
	
	Exon_ExonID=ExonLevelData[which(ExonLevelData$ExonID==ExonID),]
	Exon_ExonID=Exon_ExonID[,-c(1,2)]
	Exon_ExonID=log2(Exon_ExonID)
	Exon_ExonID=Exon_ExonID[,unlist(groups)]
	
	
	# Gene Level Data
	Gene_GeneID=GeneLevelData[which(GeneLevelData$GeneID==GeneID),]
	Gene_GeneID[,-c(1)]=log2(Gene_GeneID[,-c(1)])
	Gene_GeneID=Gene_GeneID[,-c(1)]
	Gene_GeneID=Gene_GeneID[,unlist(groups)]
	
	group1=seq_along(groups$group1)
	group2=length(groups$group1)+seq_along(groups$group2)
#	print(t.test(x=Gene_GeneID[1,group1],y=Gene_GeneID[1,group2]))
	
	
	# Observed Probe intensities
	ProbeIntensities=Data[which(Data$ExonID==ExonID),]
	ProbeIntensities_ExonID=ProbeIntensities[,-c(1,2)]
	ProbeIntensities_ExonID=ProbeIntensities_ExonID[,unlist(groups)]
	
	ProbeIntensities_ExonID=apply(ProbeIntensities_ExonID,2,as.character)
	ProbeIntensities_ExonID=apply(ProbeIntensities_ExonID,2,as.numeric)
	#plot 1 : gene level + exon level + probe intensities
	
	max1=max(Gene_GeneID,Exon_ExonID)
	max_ylim=max(as.numeric(ProbeIntensities_ExonID),max1)
	if(!is.null(location)){
		location1=paste(location,"_plot1",sep="")
	}
	else{
		location1=NULL
	}
	plottypein(plottype,location1)
	par(mar=c(6.5,4.5,2.5,4.5))
	if(is.null(ylabel)){
		ylabel=paste("Transcript",GeneID," - PSR",ExonID,sep=" ")
	}
	
	plot(0,0,typ="n",xlab="",ylab=ylabel,ylim=c(0,max_ylim+2),xlim=c(1,ncol(ProbeIntensities_ExonID)),xaxt="n")
	lines(x=c(1:ncol(ProbeIntensities_ExonID)),y=Gene_GeneID) #Gene level data...
	lines(x=c(1:ncol(ProbeIntensities_ExonID)),y=Exon_ExonID,col="blue") #Exon Level data
	for(i in 1:nrow(ProbeIntensities_ExonID)){
		points(x=c(1:ncol(ProbeIntensities_ExonID)),y=as.matrix(ProbeIntensities_ExonID[i,]),pch=19,col="blue")
	}
	axis(1,labels=colnames(Gene_GeneID),at=c(1:ncol(ProbeIntensities_ExonID)),las=2)
	plottypeout(plottype)
	
	#plot reserve : plot of cdf of exon scores
	
#	cdf=ecdf(ExonScores$X50.)
#	plot(sort(ExonScores$X50.),cdf(sort(ExonScores$X50.)),pch=19,xlim=c(0,1),ylim=c(0,1),xlab="Exon scores",ylab="Fn(x)")
#	points(sort(ExonScores$X50.[which(ExonScores$exon%in%ASExons_2736322$ExonID)]),cdf(sort(ExonScores$X50.[which(ExonScores$exon%in%ASExons_2736322$ExonID)])),pch=19,col="red")
	
	
	
}
PlotFunctionRanks<-function(GeneID=NULL,ExonID=NULL,Data,groups){
	
	if(is.null(GeneID) | is.null(ExonID)){
		stop("no GeneID and/or ExonID specified")
	}
	
	# Observed Probe intensities
	Ranks=Data[,unlist(groups)]
	
	max_ylim=max(as.vector(as.matrix(Ranks)))
	
	ylabel=paste("Transcript",GeneID," - PSR",ExonID,sep=" ")
	
	plot(0,0,typ="n",xlab="",ylab=ylabel,ylim=c(0,max_ylim+2),xlim=c(1,ncol(Ranks)),xaxt="n")
	for(i in 1:nrow(Ranks)){
		text(x=c(1:ncol(Ranks)),y=as.matrix(Ranks[i,]),label=as.character(Ranks[i,]),col="blue",cex=1.2,font=2)
	}
	axis(1,labels=colnames(Ranks),at=c(1:ncol(Ranks)),las=2)
	
}

## Genome-Wide Analysis 

length(unique(HTAData_RASA$GeneID))
# [1] 33516
length(unique(HTAData_RASA$ExonID))
# [1] 547756
AllExons=unique(HTAData_RASA$ExonID)
length(which(substr(AllExons,1,3)=="PSR"))
# [1] 298281
length(which(substr(AllExons,1,3)=="JUC"))
# [1] 249475

# Number of genes remaining after I/NI
length(unique(HTAData_RASA_PSR_LiverVSMuscle_AllPSR$geneID))
# [1] 12597

# Number of probe sets remaining after I/NI:
length(unique(HTAData_RASA_PSR_LiverVSMuscle_AllPSR$ExonID))
# [1] 148848      6

# Number of AS Probe sets:
Filter=HTAData_RASA_PSR_LiverVSMuscle_AllPSR[HTAData_RASA_PSR_LiverVSMuscle_AllPSR$X50.>=0.50,]
#Filter=HTAData_RASA_PSR_LiverVSMuscle_AllPSR[HTAData_RASA_PSR_LiverVSMuscle_AllPSR$ExonScore>=0.50,]
ASPSR_PSR=Filter[Filter$adj.p.value<0.05,]
length(unique(ASPSR_PSR$ExonID))
# [1] 20992


# Top Probe set PSR010025633
head(ASPSR_PSR)
#      geneID       ExonID      X50. t.statistic      p.value adj.p.value
# 1 TC0102569 PSR010025633 0.9865651   -30.89823 3.584204e-04 0.008546966
# 2 TC0700986 PSR070008232 0.9800008   -41.95576 1.972998e-06 0.001167280
# 3 TC1201038 PSR120009705 0.9770366   -49.37414 3.785028e-04 0.008766048
# 4 TC0301262 PSR030011880 0.9759725   -75.41089 1.189661e-04 0.005078599
# 5 TC1200691 PSR120006355 0.9744704    28.43162 3.943601e-05 0.003186866
# 6 TC0100415 PSR010004149 0.9730301   -55.69742 3.280375e-06 0.001295165

PlotFunction2(GeneID="TC0102569",ExonID="PSR010025633",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - PSR010025633",plottype="sweave",location=NULL)


# Top Probe set PSR080007308
ASPSR_PSR[which(ASPSR_PSR$ExonID=="PSR080007308"),]
#       geneID       ExonID      X50. t.statistic      p.value adj.p.value
# 68 TC0800969 PSR080007308 0.9432664   -30.81983 0.0006971904  0.01161996


PlotFunction2(GeneID="TC0800969",ExonID="PSR080007308",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)



# Distribution of the 20 992 identified exons

load("Data/REIDS_RankedProbesets.RData")
REIDS_RankedProbesets=data.frame(REIDS_RankedProbesets)
for(i in c(1:6)){
	REIDS_RankedProbesets[,i]=as.character(REIDS_RankedProbesets[,i])
}
REIDS_RankedProbesets[,4]=as.numeric(as.character(REIDS_RankedProbesets[,4]))
REIDS_RankedProbesets[,5]=as.numeric(as.character(REIDS_RankedProbesets[,5]))
REIDS_RankedProbesets=REIDS_RankedProbesets[order(REIDS_RankedProbesets[,5]),]
head(REIDS_RankedProbesets)
#              X1            X2            X3 X4 X5                      X6
# 1  PSR010025633 JUC0100236357 JUC0100236756  1  0                        
# 6  PSR010004149 JUC0100031528 JUC0100034464  1  0                        
# 7  PSR170010524 JUC1700094642                1  0 Only 1 type of junction
# 10 PSR020002230 JUC0200020846                1  0 Only 1 type of junction
# 17 PSR010007541 JUC0100065769                1  0 Only 1 type of junction
# 36 PSR010013561 JUC0100113569                1  0 Only 1 type of junction


dim(REIDS_RankedProbesets)
# [1] 20992     6
## probe sets with 2 linking junctions:
length(which(REIDS_RankedProbesets[,6]==""|REIDS_RankedProbesets[,6]=="1 of two junctions is supporting"))
# [1] 10838
## probe sets with 2 linking junctions and both are supporting
length(which(REIDS_RankedProbesets[,6]==""&round(REIDS_RankedProbesets[,4],2)>=0.05))
# [1] 6658
## probe sets with 2 linking junctions but has a sign interaction term
length(which((REIDS_RankedProbesets[,6]==""|REIDS_RankedProbesets[,6]=="1 of two junctions is supporting")&round(REIDS_RankedProbesets[,4],2)<0.05))
# [1] 4180
## probe sets with 2 linking junctions and 1 is supporting
length(which(REIDS_RankedProbesets[,6]=="1 of two junctions is supporting"&round(REIDS_RankedProbesets[,4],2)<0.05))
# [1] 2888
## probe sets with 2 linking junctions and 1 is supporting
length(which(REIDS_RankedProbesets[,6]==""&round(REIDS_RankedProbesets[,4],2)<0.05))
# [1] 1292
## probe sets with 1 or no linking junctions
length(which(REIDS_RankedProbesets[,6]!=""&REIDS_RankedProbesets[,6]!="1 of two junctions is supporting"))
# [1] 10154
## probe sets without linking junctions
length(which(REIDS_RankedProbesets[,6]=="Only exclusion junctions"))
# [1] 2981
## probe sets without linking junctions
length(which(REIDS_RankedProbesets[,6]=="No Annot Junction"))
# [1] 2110
## probe sets with 1 linking junction
length(which(REIDS_RankedProbesets[,6]=="Only 1 type of junction"))
# [1] 5062
## probe sets with 1 linking junction and supporting
length(which(REIDS_RankedProbesets[,6]=="Only 1 type of junction"&REIDS_RankedProbesets[,4]>0.05))
# [1] 3455
## probe sets with 1 linking junction but not supporting
length(which(REIDS_RankedProbesets[,6]=="Only 1 type of junction"&REIDS_RankedProbesets[,4]<0.05))
# [1] 1607

save(REIDS_RankedProbesets,file="Data/REIDS_RankedProbesets.RData")



## Alternative Splicing With Supporting Junctions

ASPSR_PSR[ASPSR_PSR$ExonID=="PSR010025633",]
#      geneID       ExonID      X50. t.statistic      p.value adj.p.value
# 1 TC0102569 PSR010025633 0.9865651   -30.89823 0.0003584204 0.008546966

m <- matrix(c(1,1,2,2,3,3,4,5,5,6,6,7),nrow = 2,ncol = 6,byrow = TRUE)
layout(mat = m)
#m <- matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
#layout(mat = m,heights=c(0.08,0.92))
#par(mar=c(0,0,0,0))
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#legend(x = "top",inset=0,legend = c("Observed Probe Values", "Probeset Level", "Transcript Level"), col=c("blue","blue","black"), lty=c(NA,1,2),pch=c(19,NA,NA), lwd=2, cex=2.2, horiz = TRUE,bty = "n",xjust=0.5)
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC0102569",ExonID="JUC0100236756",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - JUC0100236756",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102569",ExonID="PSR010025633",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - PSR010025633",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102569",ExonID="JUC0100236357",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - JUC0100236357",plottype="sweave",location=NULL)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC0102569",ExonID="PSR010025632",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569- PSR010025632",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102569",ExonID="PSR010025634",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - PSR010025634",plottype="sweave",location=NULL)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

par(mfrow=c(1,3))
PlotFunction2(GeneID="TC0102569",ExonID="JUC0100236286",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569- JUC0100236286",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102569",ExonID="JUC0100236573",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - JUC0100236573",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102569",ExonID="JUC0100236628",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - JUC0100236628",plottype="sweave",location=NULL)



#Figures with observed values in upper panel and conversion to ranks in lower panel
#Conversion to ranks (taken from function)
#Supplmentary Figure 4
Data=HTAData_RASA
ASPSR="PSR010025633"
Juc3="JUC0100236756"
Juc5="JUC0100236357"
PSR=as.vector(as.matrix(Data[which(Data$ExonID==as.character(ASPSR)),-c(1,2)]))
PSR_Ranks=sort(PSR,index.return=TRUE)$ix
JUC3=as.vector(as.matrix(Data[which(Data$ExonID==Juc3),-c(1,2)]))
JUC3_Ranks=sort(JUC3,index.return=TRUE)$ix
JUC5=as.vector(as.matrix(Data[which(Data$ExonID==Juc5),-c(1,2)]))
JUC5_Ranks=sort(JUC5,index.return=TRUE)$ix

J5=matrix(JUC5_Ranks,ncol=6,nrow=8)
colnames(J5)=c("Liver_A","Liver_B","Liver_C","Muscle_A","Muscle_B","Muscle_C")
J3=matrix(JUC3_Ranks,ncol=6,nrow=8)
colnames(J3)=c("Liver_A","Liver_B","Liver_C","Muscle_A","Muscle_B","Muscle_C")
P=matrix(PSR_Ranks,ncol=6,nrow=8)
colnames(P)=c("Liver_A","Liver_B","Liver_C","Muscle_A","Muscle_B","Muscle_C")

L=list(PSR_Ranks,JUC3_Ranks,JUC5_Ranks)
if(length(unique(c(length(PSR),length(JUC3),length(JUC5))))>1){
	MinLength=min(c(length(PSR),length(JUC3),length(JUC5)))
	Index=which(c(length(PSR),length(JUC3),length(JUC5))>MinLength)
	for(s in Index){
		L[[s]]=L[[s]][-which(L[[s]]%in%c((MinLength+1):length(L[[s]])))]
	}
}

Ranks=cbind(L[[3]],L[[1]],L[[2]])

Y=as.vector(as.matrix(Ranks))
exon=c(rep(1,nrow(Ranks)),rep(2,nrow(Ranks)),rep(3,nrow(Ranks)))
tissue=rep(c(rep(1,nrow(Ranks)/2),rep(2,nrow(Ranks)/2)),3)

exon=as.factor(exon)
tissue=as.factor(tissue)
ft1<-lm(Y~exon*tissue-1)
ft2<-lm(Y~exon+tissue)

Chisq=lrtest(ft2,ft1)[2,4]


m <- matrix(c(1,1,2,2,3,3,4,4,5,5,6,6),nrow = 2,ncol = 6,byrow = TRUE)
layout(mat = m)
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC0102569",ExonID="JUC0100236756",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - JUC0100236756",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102569",ExonID="PSR010025633",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - PSR010025633",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102569",ExonID="JUC0100236357",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0102569 - JUC0100236357",plottype="sweave",location=NULL)
PlotFunctionRanks(GeneID="TC0102569",ExonID="JUC0100236756",Data=J3,groups=list(group1=c(1,2,3),group2=c(4,5,6)))
PlotFunctionRanks(GeneID="TC0102569",ExonID="PSR010025633",Data=P,groups=list(group1=c(1,2,3),group2=c(4,5,6)))
PlotFunctionRanks(GeneID="TC0102569",ExonID="JUC0100236756",Data=J5,groups=list(group1=c(1,2,3),group2=c(4,5,6)))



## Alternative Splicing Without Supporting Junctions

which(REIDS_RankedProbesets[,1]=="PSR080007308")
# [1] 15787
REIDS_RankedProbesets[REIDS_RankedProbesets$X1=="PSR080007308",]
#              X1            X2            X3           X4       X5 X6
# 67 PSR080007308 JUC0800063802 JUC0800064097 3.193896e-10 46.85273   

#5' and 3'
par(mfrow=c(1,4))
PlotFunction2(GeneID="TC0800969",ExonID="JUC0800064097",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969 - JUC0800064097",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0800969",ExonID="PSR080007308",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969 - PSR080007308",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0800969",ExonID="JUC0800063802",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969  - JUC0800063802",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0800969",ExonID="PSR080007307",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969  - PSR080007307",plottype="sweave",location=NULL)
#exclusion
par(mfrow=c(1,3))
PlotFunction2(GeneID="TC0800969",ExonID="JUC0800063922",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969 - JUC0800063922",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0800969",ExonID="JUC0800064078",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969 - JUC0800064078",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0800969",ExonID="JUC0800064160",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969 - JUC0800064160",plottype="sweave",location=NULL)

#PSR080007309 and PSR080007310
par(mfrow=c(1,2))
PlotFunction2(GeneID="TC0800969",ExonID="PSR080007309",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969 - PSR080007309",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0800969",ExonID="PSR080007310",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC0800969 - PSR080007310",plottype="sweave",location=NULL)


## Alternative Splicing with a "single side" non Supporting Junction

# PSR140003000
HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$PSR_ID=="PSR140003000"),]
#            TC_ID       PSR_ID  PSR_X50. PSR_t.statistic  PSR_p.value
# 226326 TC1400365 PSR140003000 0.9538094        -73.1473 4.138308e-07
#        PSR_adj.p.value     psr_type        JUC_ID  JUC_X50. JUC_t.statistic
# 226326    0.0003903692 constitutive JUC1400033138 0.6594335        7.229892
#        JUC_p.value JUC_adj.p.value junction_type junction_probeset_id as_type
# 226326 0.002372053      0.01615952  constitutive               354667       3
#             Status
# 226326 misbehaving

par(mfrow=c(2,3))
PlotFunction2(GeneID="TC1400365",ExonID="PSR140003000",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1400365",ExonID="JUC1400033138",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1400365",ExonID="PSR140003001",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)


which(REIDS_RankedProbesets[,1]=="PSR140003000")
# [1] 15765



## RASA and REIDS

load(file="Data/RASA_ASPSR.RData")
Mapping=read.table("F:/Doctorate/Detection of Splicing Variants/Preprocessing/Aroma Package/HTA 6 samples - run 2/Creating CDF/TC_PSR_JUCID_NumberID.txt",sep="\t",header=FALSE)
RASAID=RASA_ASPSR$asa
colnames(Mapping)=c("PSR","PSR_ID","TC_extra")
RASA=merge(RASAID,Mapping,by="PSR")

length(unique(RASA$PSR_ID))
# [1] 49014
length(unique(REIDS_RankedProbesets[,1]))
# [1] 20992
length(intersect(REIDS_RankedProbesets[,1],RASA$PSR_ID))
# [1] 9196 
length(intersect(REIDS_RankedProbesets[which(REIDS_RankedProbesets[,6]==""&round(REIDS_RankedProbesets[,4]>=0.05,2)),1],RASA$PSR_ID))
# [1] 3444
length(intersect(REIDS_RankedProbesets[which(REIDS_RankedProbesets[,4]!=100000),1],RASA$PSR_ID))
# [1] 7763


# ranks
par(mar=c(4.5,4.5,0.3,0.2))
hist(which(REIDS_RankedProbesets[which(REIDS_RankedProbesets[,4]!=100000),1]%in%RASA$PSR_ID),main="",xlab="Ranks",n=20)


#PSR140003960
m <- matrix(c(1,1,2,2,3,3,4,4,5,5,6,6),nrow = 2,ncol = 6,byrow = TRUE)
layout(mat = m)
#m <- matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
#layout(mat = m,heights=c(0.08,0.92))
#par(mar=c(0,0,0,0))
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#legend(x = "top",inset=0,legend = c("Observed Probe Values", "Probeset Level", "Transcript Level"), col=c("blue","blue","black"), lty=c(NA,1,2),pch=c(19,NA,NA), lwd=2, cex=2.2, horiz = TRUE,bty = "n",xjust=0.5)
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC1400481",ExonID="PSR140003960",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1400481",ExonID="JUC1400040768",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1400481",ExonID="JUC1400040831",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)

#Identified by RASA (not by REIDS)
unique(RASA$PSR_ID[which(!RASA$PSR_ID%in%REIDS_RankedProbesets[,1])])[c(1:5)]

HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$PSR_ID=="PSR030013789"),c(1,2,8,15,16)]
#           TC_ID       PSR_ID        JUC_ID as_type              Status
# 73015 TC0301476 PSR030013789 JUC0300131708       5                  ok
# 73019 TC0301476 PSR030013789 JUC0300131720       3         misbehaving
# 73026 TC0301476 PSR030013789 JUC0300131735       5 filtered out by INI
# 73040 TC0301476 PSR030013789 JUC0300131773       3 filtered out by INI
# 73042 TC0301476 PSR030013789 JUC0300131776       5         misbehaving

HTAData_RASA_PSR_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSR_LiverVSMuscle_AllPSR$ExonID=="PSR030013789"),]
#           geneID       ExonID      X50. t.statistic      p.value adj.p.value
# 127245 TC0301476 PSR030013789 0.3456932   -26.03578 3.020016e-05 0.002860025

m <- matrix(c(1,1,2,2,3,3,4,4,5,5,6,6),nrow = 2,ncol = 6,byrow = TRUE)
layout(mat = m)
#m <- matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
#layout(mat = m,heights=c(0.08,0.92))
#par(mar=c(0,0,0,0))
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#legend(x = "top",inset=0,legend = c("Observed Probe Values", "Probeset Level", "Transcript Level"), col=c("blue","blue","black"), lty=c(NA,1,2),pch=c(19,NA,NA), lwd=2, cex=2.2, horiz = TRUE,bty = "n",xjust=0.5)
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC0301476",ExonID="JUC0300131720",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0301476",ExonID="PSR030013789",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0301476",ExonID="JUC0300131708",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)


#Identified by REIDS (not by RASA)
REIDS_RankedProbesets[19,]
#               X1            X2            X3 X4 X5
# 307 PSR050011925 JUC0500105779 JUC0500105725  0   

ASPSR_PSRwithJUN[ASPSR_PSRwithJUN$PSR_ID=="PSR050011925",c(1,2,8,15,16)]
#            TC_ID       PSR_ID        JUC_ID as_type Status
# 101182 TC0501406 PSR050011925 JUC0500105725       5     ok
# 101192 TC0501406 PSR050011925 JUC0500105779       3     ok

m <- matrix(c(1,1,2,2,3,3,4,5,6,6,7,8),nrow = 2,ncol = 6,byrow = TRUE)
layout(mat = m)
#m <- matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
#layout(mat = m,heights=c(0.08,0.92))
#par(mar=c(0,0,0,0))
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#legend(x = "top",inset=0,legend = c("Observed Probe Values", "Probeset Level", "Transcript Level"), col=c("blue","blue","black"), lty=c(NA,1,2),pch=c(19,NA,NA), lwd=2, cex=2.2, horiz = TRUE,bty = "n",xjust=0.5)
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC0501406",ExonID="JUC0500105725",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0501406",ExonID="PSR050011925",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0501406",ExonID="JUC0500105779",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

HTAData_RASA_PSR_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSR_LiverVSMuscle_AllPSR$ExonID=="PSR050011925"),]
#        geneID       ExonID      X50. t.statistic      p.value adj.p.value
# 313 TC0501406 PSR050011925 0.9004739   -34.97897 0.0001628763  0.00590736


# Supported by one junction
REIDS_RankedProbesets[which(REIDS_RankedProbesets[,1]=="PSR010025473"),]
#                X1            X2            X3         X4       X5
# 5913 PSR010025473 JUC0100235426 JUC0100235442 0.04145239 6.643221
#                                    X6
# 5913 1 of two junctions is supporting

HTAData_RASA_PSR_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSR_LiverVSMuscle_AllPSR$ExonID=="PSR010025473"),]
#         geneID       ExonID     X50. t.statistic     p.value adj.p.value
# 9041 TC0102559 PSR010025473 0.690093    6.121972 0.003670564  0.02654404
HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$PSR_ID=="PSR030012131"),c(1,2,8,15,16)]
#           TC_ID       PSR_ID        JUC_ID as_type      Status
# 71304 TC0301284 PSR030012131 JUC0300116077       5          ok
# 71312 TC0301284 PSR030012131 JUC0300116090       3 misbehaving
par(mfrow=c(2,3))
PlotFunction2(GeneID="TC0102559",ExonID="JUC0100235426",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102559",ExonID="PSR010025473",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0102559",ExonID="JUC0100235442",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)



#verified probesets
load("Data/VerifiedProbes_PSRID.RData")
nrow(VerifiedProbes_PSRID)
# [1] 373

length(which(as.character(VerifiedProbes_PSRID[,2])%in%unique(as.character(REIDS_RankedProbesets[,1]))))
# [1] 294
length(which(as.character(VerifiedProbes_PSRID[,2])%in%unique(as.character(REIDS_RankedProbesets[,1]))))
# [1] 294

length(which(as.character(VerifiedProbes_PSRID[,2])%in%unique(as.character(REIDS_RankedProbesets[which(REIDS_RankedProbesets[,4]!=100000),][,1]))))
# [1] 147
length(which(as.character(VerifiedProbes_PSRID[,2])%in%unique(as.character(REIDS_RankedProbesets[which(REIDS_RankedProbesets[,4]!=100000),][,1]))))
# [1] 248


Categorize<-function(Find,Data=HTAData_RASA_PSR_LiverVSMuscle_AllPSR,GDSOutput=HTAData_RASA_OutputGDSModel,SignificanceLevel=0.05,ICCThreshold=0.50,SuppJuns=TRUE,MisbJuns=FALSE){
	
	Out=c()
	
	Subset=Data[which(as.character(Data$PSR_ID)%in%as.character(Find)),]
	RASA_PSR_NotPresent=as.character(Find)[which(!as.character(Find)%in%as.character(unique(Subset$PSR_ID)))]
	
	
	for(i in 1:length(GDSOutput)){
		print(i)
		Temp=GDSOutput[[i]]$Informative
		if(length(which(as.character(RASA_PSR_NotPresent)%in%as.character(Temp$exonNames)))>0){
			TempPSR=as.character(RASA_PSR_NotPresent)[which(as.character(RASA_PSR_NotPresent)%in%as.character(Temp$exonNames))]
			for(j in 1:length(TempPSR)){
				if(Temp[as.character(Temp$exonNames)==as.character(TempPSR[j]),3]==FALSE){
					Out=rbind(Out,c(as.character(TempPSR[j]),"filtered out by INI"))
					RASA_PSR_NotPresent=as.character(RASA_PSR_NotPresent)[-which(as.character(RASA_PSR_NotPresent)==as.character(TempPSR[j]))]
				}
			}
		}
	}
	
	if(length(RASA_PSR_NotPresent)>0){
		Out=rbind(Out,cbind(as.character(RASA_PSR_NotPresent),"PSR not analyzed by GDS"))
	}
	
	colnames(Out)=c("PSR","GDS Status")
	
	RASA_PSR_Present=unique(Subset$PSR_ID)
	
	for(i in 1:length(RASA_PSR_Present)){
		print(i)
		TempGDS=Subset[which(as.character(Subset$PSR_ID)==as.character(RASA_PSR_Present[i])),]
		
		#check the different caveats
		
		# 1) ICC of PSR
		if(TempGDS$PSR_X50.[1]<ICCThreshold){
			Out=rbind(Out,c(as.character(RASA_PSR_Present[i]),"PSR does not pass ICC"))
		}
		# 2) p-value of PSR
		else if(TempGDS$PSR_adj.p.value[1]>SignificanceLevel){
			Out=rbind(Out,c(as.character(RASA_PSR_Present[i]),"PSR does not pass significancelevel"))
		}
		
		else if(SuppJuns==TRUE){
			
			Filter0=TempGDS[!is.na(TempGDS$JUC_X50.),]
			if(nrow(Filter0)==0){
				Out=rbind(Out,c(RASA_PSR_Present[i],"No junctions pass INI"))
			}
			else{
				Filter1=Filter0
				if(nrow(Filter1)==0){
					Out=rbind(Out,c(RASA_PSR_Present[i],"No junctions pass ICC"))
				}
				
				else{
					Filter2=Filter1[Filter1$JUC_adj.p.value<=SignificanceLevel,]
					if(nrow(Filter2)==0){
						Out=rbind(Out,c(RASA_PSR_Present[i],"No junctions pass significancelevel"))
					}
					else{
						if(MisbJuns==TRUE){
							Filter3=Filter2[Filter2$Status=="ok",]
							if(nrow(Filter3)==0){
								Out=rbind(Out,c(RASA_PSR_Present[i],"Misbehaving junctions"))
							}
							else{
								Out=rbind(Out,c(RASA_PSR_Present[i],"PSR identified"))
							}
						}
						else{
							Out=rbind(Out,c(RASA_PSR_Present[i],"PSR identified"))
						}
					}
				}
			}
		}
		else{
			Out=rbind(Out,c(as.character(RASA_PSR_Present[i]),"PSR identified"))
		}
	}
	
	if(nrow(Out)!=length(Find)){
		NotFound=Find[which(!Find%in%Out[,1])]
		message("returned object does not have the same length as Find.")
		return(NotFound)
	}
	Out=as.data.frame(Out)
	NumberPSRNotPassINI=length(which(Out$"GDS Status"=="filtered out by INI"))
	NumberPSRNotAnalyzed=length(which(Out$"GDS Status"=="PSR not analyzed by GDS"))
	NumberPSRNotPassICC=length(which(Out$"GDS Status"=="PSR does not pass ICC"))
	NumberPSRNotPassSign=length(which(Out$"GDS Status"=="PSR does not pass significancelevel"))
	NumberPSRNoJunsINI=length(which(Out$"GDS Status"=="No junctions pass INI"))
	NumberPSRNoJunsICC=length(which(Out$"GDS Status"=="No junctions pass ICC"))
	NumberPSRNoJunsSign=length(which(Out$"GDS Status"=="No junctions pass significancelevel"))
	NumberPSRNoJunsMisb=length(which(Out$"GDS Status"=="Misbehaving junctions"))
	NumberPSRIdentified=length(which(Out$"GDS Status"=="PSR identified"))
	
	cat(paste("Number of Probes not passing INI :",NumberPSRNotPassINI,"\n",sep=" "))
	cat(paste("Number of Probes not analyzed by GDS :",NumberPSRNotAnalyzed,"\n",sep=" "))
	cat(paste("Number of Probes not passing ICC :",NumberPSRNotPassICC,"\n",sep=" "))
	cat(paste("Number of Probes not passing significancelevel :",NumberPSRNotPassSign,"\n",sep=" "))
	cat(paste("Number of Probes without juns passing INI :",NumberPSRNoJunsINI,"\n",sep=" "))
	cat(paste("Number of Probes without juns passing ICC :",NumberPSRNoJunsICC,"\n",sep=" "))
	cat(paste("Number of Probes without juns passing significancelevel :",NumberPSRNoJunsSign,"\n",sep=" "))
	cat(paste("Number of Probes with misbehaving junctions :",NumberPSRNoJunsMisb,"\n",sep=" "))
	cat(paste("Number of identified Probes :",NumberPSRIdentified,"\n",sep=" "))
	
	return(Out)
	
	
	
}

FindRasaPSR=Categorize(Find=as.character(VerifiedProbes_PSRID$id))

which(REIDS_RankedProbesets[,1]%in%VerifiedProbes_PSRID[,2])

t=HTAData_RASA_PSR_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSR_LiverVSMuscle_AllPSR$ExonID%in%VerifiedProbes_PSRID[,2]),]
length(which(t[,3]<0.50))
# [1] 43
length(which(t[,3]>0.50&t[,6]<0.05))
# [1] 294
length(which(t[,3]>0.50&t[,6]>0.05))
# [1] 13
t[which(t[,3]>0.50&t[,6]>0.05),]
#          geneID       ExonID      X50. t.statistic    p.value adj.p.value
# 3613  TC0301825 PSR030016748 0.7652603   -5.109297 0.02538637  0.08061593
# 5221  TC0500702 PSR050006184 0.7373206    6.024853 0.02232692  0.07434480
# 8144  TC1900984 PSR190009115 0.6995667   -7.144133 0.01446133  0.05730954
# 11644 TC0801225 PSR080009043 0.6688895    7.122764 0.01290099  0.05363479
# 13218 TC0500530 PSR050004757 0.6578463   -6.852927 0.01923468  0.06792349
# 28288 TC0202331 PSR020021147 0.5870701    4.555179 0.01215476  0.05173464
# 30365 TC0202295 PSR020020801 0.5800056   -5.642821 0.01867651  0.06672175
# 40143 TC0200358 PSR020003217 0.5520871    4.426533 0.04647691  0.11855424
# 47182 TC0900703 PSR090006386 0.5347544   -7.404050 0.01509913  0.05878054
# 49276 TC1400248 PSR140001763 0.5299351   -4.895132 0.01283155  0.05344463
# 50826 TC1600350 PSR160003823 0.5264922   -5.567754 0.02087982  0.07130872
# 52432 TC1200068 PSR120000756 0.5229954    3.198140 0.06150580  0.14250756
# 60732 TC0100199 PSR010002050 0.5063221    3.359816 0.04190765  0.11077336



## Supplementary Material 

#Junction Design


#PSR190012156

par(mfrow=c(2,2))
PlotFunction2(GeneID="TC1901243",ExonID="JUC1900096693",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - JUC1900096693",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1901243",ExonID="JUC1900096711",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - JUC1900096711",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1901243",ExonID="PSR190012156",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012156",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1901243",ExonID="PSR190012157",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012157",plottype="sweave",location=NULL)


m <- matrix(c(1,2,2,3,3,4,4,5,6,6,7,7,8,8,9,9),nrow = 2,ncol = 8,byrow = TRUE)
layout(mat = m)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC1901243",ExonID="PSR190012156",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012156",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1901243",ExonID="JUC1900096716",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - JUC1900096711",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1901243",ExonID="PSR190012150",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012150",plottype="sweave",location=NULL)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC1901243",ExonID="PSR190012155",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012155",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1901243",ExonID="PSR190012153",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012153",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1901243",ExonID="PSR190012152",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012152",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1901243",ExonID="PSR190012151",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012151",plottype="sweave",location=NULL)




# PSR120012240

HTAData_RASA_PSR_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSR_LiverVSMuscle_AllPSR$ExonID=="PSR120012240"),]

HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$PSR_ID=="PSR120012240"),]

HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$JUC_ID=="JUC1200112452"),]

par(mfrow=c(1,3))
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012240",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012240",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1201288",ExonID="JUC1200112452",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - JUC1200112452",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012241",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012241",plottype="sweave",location=NULL)

HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$JUC_ID=="JUC1200112431"),]

m <- matrix(c(1,1,2,2,3,3,4,5,5,6,6,7),nrow = 2,ncol = 6,byrow = TRUE)
layout(mat = m)
#m <- matrix(c(1,1,2,3),nrow = 2,ncol = 2,byrow = TRUE)
#layout(mat = m,heights=c(0.08,0.92))
#par(mar=c(0,0,0,0))
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#legend(x = "top",inset=0,legend = c("Observed Probe Values", "Probeset Level", "Transcript Level"), col=c("blue","blue","black"), lty=c(NA,1,2),pch=c(19,NA,NA), lwd=2, cex=2.2, horiz = TRUE,bty = "n",xjust=0.5)
par(mar=c(6.5,4.5,2.5,4.5))
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012237",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012237",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1201288",ExonID="JUC1200112431",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - JUC1200112431",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012240",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012240",plottype="sweave",location=NULL)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012238",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012238",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012239",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012239",plottype="sweave",location=NULL)
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$JUC_ID=="JUC1200112702"),]

par(mfrow=c(2,2))
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012238",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012238",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1201288",ExonID="JUC1200112702",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - JUC1200112702",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012240",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012240",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1201288",ExonID="PSR120012239",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1201288 - PSR120012239",plottype="sweave",location=NULL)

## Examples of Probe Sets Supported by All Linking Junctions

#PSR030011882
HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$PSR_ID=="PSR030011882"),]

par(mfrow=c(2,3))
PlotFunction2(GeneID="TC0301262",ExonID="JUC0300113792",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - JUC1900096711",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0301262",ExonID="PSR030011882",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - PSR190012156",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0301262",ExonID="JUC0300114121",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - JUC1900096716",plottype="sweave",location=NULL)

PlotFunction2(GeneID="TC0301262",ExonID="JUC0300113827",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - JUC1900096716",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0301262",ExonID="PSR030011889",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - JUC1900096716",plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC0301262",ExonID="JUC0300113838",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel="TC1901243 - JUC1900096716",plottype="sweave",location=NULL)



#TC1601187
load("Data/VerifiedProbes_PSRID.RData")
VerifiedProbes_PSRID[which(VerifiedProbes_PSRID$transcript_cluster_id=="TC1601187"),]
#         probeset_id           id transcript_cluster_id
# 728126       442409 PSR160012163             TC1601187
# 754877       458741 PSR160012166             TC1601187
# 913753       555336 PSR160012168             TC1601187
# 1073725      652605 PSR160012170             TC1601187
# 1090229      662580 PSR160012165             TC1601187

HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$PSR_ID%in%VerifiedProbes_PSRID[which(VerifiedProbes_PSRID$transcript_cluster_id=="TC1601187"),2]),]



m <- matrix(c(1,1,2,2,3,3,4,5,5,6,6,7),nrow = 2,ncol = 6,byrow = TRUE)
layout(mat = m)
par(mar=c(6.5,4.5,2.5,4.5))
for(i in VerifiedProbes_PSRID[which(VerifiedProbes_PSRID$transcript_cluster_id=="TC1601187"),2][c(1:3)]){
	PlotFunction2(GeneID="TC1601187",ExonID=i,Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel=NULL,plottype="sweave",location=NULL)
}
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
par(mar=c(6.5,4.5,2.5,4.5))
for(i in VerifiedProbes_PSRID[which(VerifiedProbes_PSRID$transcript_cluster_id=="TC1601187"),2][c(4,5)]){
	PlotFunction2(GeneID="TC1601187",ExonID=i,Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel=NULL,plottype="sweave",location=NULL)
}
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")


HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$PSR_ID=="PSR160012167"),]

par(mfrow=c(2,3))
PlotFunction2(GeneID="TC1601187",ExonID="JUC1600099690",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1601187",ExonID="PSR160012167",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1601187",ExonID="JUC1600100109",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction3(GeneID="TC1601187",ExonID="JUC1600099690",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction3(GeneID="TC1601187",ExonID="PSR160012167",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction3(GeneID="TC1601187",ExonID="JUC1600100109",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)



HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR[which(HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR$PSR_ID=="PSR160012181"),]

par(mfrow=c(2,4))
PlotFunction2(GeneID="TC1601187",ExonID="JUC1600100261",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel=paste("TC1601187 - The The 5", expression("\'"), "Junction JUC1600100261",sep=""),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1601187",ExonID="PSR160012181",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel=paste("TC1601187 - The PSR PSR160012181",sep=""),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1601187",ExonID="JUC1600099680",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel=paste("TC1601187 - The The 5", expression("\'"), "Junction JUC1600099680",sep=""),plottype="sweave",location=NULL)
PlotFunction2(GeneID="TC1601187",ExonID="JUC1600100154",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),ylabel=paste("TC1601187 - The The 5", expression("\'"), "Junction JUC1600100154",sep=""),plottype="sweave",location=NULL)

PlotFunction3(GeneID="TC1601187",ExonID="JUC1600100261",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction3(GeneID="TC1601187",ExonID="PSR160012181",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction3(GeneID="TC1601187",ExonID="JUC1600099680",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)
PlotFunction3(GeneID="TC1601187",ExonID="JUC1600100154",Data=HTAData_RASA,GDS_Output=HTAData_RASA_OutputGDSModel,GeneLevelData=HTA_GeneLevelSummarized_rma_RASA,ExonLevelData=HTA_ExonLevelSummarized_rma_RASA,FIRMA_Data=NULL,groups=list(group1=c(1,2,3),group2=c(4,5,6)),plottype="sweave",location=NULL)

