
#Data
load("HTAData_RASA_WithoutJunctions_OutputREIDSModel.RData")
Probes_and_Junctions=read.table(file="hta1.r1.ass",sep="\t",header=TRUE)

#Groups
groupLiver=c(1,2,3)
groupMuscle=c(4,5,6)

ListofGroups=list(groupLiver,groupMuscle)
names(ListofGroups)=c('groupLiver','groupMuscle')

#Function

# Perfoms t-test on the array scores and provides a data frame with junction information  
ExonTesting <- function(Data, Exonthreshold=NULL,groups=list(group1=NULL,group2=NULL),paired=FALSE,significancelevel=NULL){
	
	if(is.null(Exonthreshold)){
		Exonthreshold=0
	}
	
	if(is.null(groups$group1)){
		message("A test on the array scores will NOT be performed")
	}
	
	message("Data is filtered to only contain genes with informative exons")
	Data_Filtered1=FilterInformativeGenes(Data)  #only those genes with informative exons remain
	
	
	filterASExons<-function(i,geneID,Subset,Exon,groupings,pairing){
		
		## Step 1 : filtering on the Exon value. If 0, no filtering occurs
		ExonsExon = Subset$exonScore
		ArrayScores = Subset$arrayScore
		
		SelectExonsExon = ExonsExon[which(ExonsExon$X50>Exon),]
		ExonsPassedExon = SelectExonsExon$exon
		
		SelectRowsArray = which(rownames(ArrayScores)%in%ExonsPassedExon)
		
		## Step 2 : Testing of the Array Scores -- paired or not paired
		if(pairing==FALSE){  # Test between two groups of Array Scores
			
			if(!(is.null(groupings$group1)) & length(SelectRowsArray)!=0){
				
				ArrayScore_group1=ArrayScores[SelectRowsArray,groupings$group1,drop=FALSE]
				ArrayScore_group2=ArrayScores[SelectRowsArray,groupings$group2,drop=FALSE]
				
				ttest<-function(i,g1,g2,pairs){
					out1=t.test(x=g1,y=g2)
					out2=cbind(out1$statistic,out1$p.value)
					
					return(out2)
				}
				
				ArrayScoreTest = t(sapply(c(1:length(ExonsPassedExon)),function(i) ttest(g1=ArrayScore_group1[i,],g2=ArrayScore_group2[i,], pairs = pairing) ))
				ArrayScoreTest=as.data.frame(ArrayScoreTest)
				colnames(ArrayScoreTest)=c("t.statistic","p.value")
				rownames(ArrayScoreTest)=rownames(ArrayScore_group1)
				
				
			}
			else{
				ArrayScoreTest = NULL
			}
		}
		
		
		if(pairing==TRUE){ # Test the mean paired difference against zero
			
			
			if(!(is.null(groupings$group1)) & length(SelectRowsArray)!=0){
				
				ArrayScore_group1=ArrayScores[SelectRowsArray,groupings$group1,drop=FALSE]
				ArrayScore_group2=ArrayScores[SelectRowsArray,groupings$group2,drop=FALSE]
				
				mean_paired_diff<-function(g1,g2){
					Paired_Diff=g1-g2
					Mean_Diff=mean(Paired_Diff)
					
					out1=t.test(x=Paired_Diff)
					out2=cbind(out1$statistic,out1$p.value)
					out3=cbind(Mean_Diff,out2)
					
					return(out3)
				}
				
				ArrayScoreTest = t(sapply(c(1:length(ExonsPassedExon)),function(i) mean_paired_diff(g1=ArrayScore_group1[i,],g2=ArrayScore_group2[i,]) ))
				ArrayScoreTest=data.frame("Mean_Diff"=ArrayScoreTest[,1],"t-statistic"=ArrayScoreTest[,2],p.value=ArrayScoreTest[,3])
				rownames(ArrayScoreTest)=rownames(ArrayScore_group1)
				
				
			}
			else{
				ArrayScoreTest = NULL
			}
		}	
		
		if(!is.null(ArrayScoreTest)){
			Output=cbind("geneID"=rep(geneID,length(ExonsPassedExon)), SelectExonsExon,ArrayScoreTest)			
		}
		else{
			Output=cbind("geneID"=rep(geneID,length(ExonsPassedExon)),SelectExonsExon)
		}
		
		colnames(Output)[2]="ExonID"
		
		if(nrow(Output)==0){
			Output=NULL
		}
		
		return(Output)
		
	}
	
	Out1=lapply(1:length(Data_Filtered1), function(i) filterASExons(i,geneID = names(Data_Filtered1[i]), Subset = Data_Filtered1[[i]], Exon = Exonthreshold , groupings=groups, pairing = paired))
	names(Out1)=names(Data_Filtered1)
	
	Out2=Out1[vapply(Out1, Negate(is.null), NA)]
	
	Out3<-do.call(rbind.data.frame, Out2)
	
	
	
	## Step 3 : if Exon threshold specified: we adjust for multiplicity ;  filter on significance level 
	if(nrow(Out3)!=0){
		message("The p-value are adjusted")
		if(!is.null(Out3$p.value)){
			Out3$adj.p.value=p.adjust(Out3$p.value,"fdr")
		}	
		if(!is.null(significancelevel)){
			Out3=Out3[which(Out3$adj.p.value<=significancelevel),]
		}
		rownames(Out3)=seq(1:nrow(Out3))
	}
	
	return(Out3)
	
}
FilterInformativeGenes <- function(Data){
	Out1=list()
	filter <- function(Subset){
		if("exonScore"%in% names(Subset)){
			return(Subset)
		}
	}
	Out1=lapply(Data,filter)
	
	Out2=Out1[vapply(Out1, Negate(is.null), NA)]
	
	return(Out2)
}
AS_PSR_JUN<-function(Data,Exonthreshold=0.5,groups=list(group1=NULL,group2=NULL),paired=FALSE,significancelevel=0.05,JunctionInfoFile=NULL){
	
	if(class(Data)=="list"){
		message("The data is assumed to be output of the GDS model. Testing the array scores will be performed")
		#message("The used threshold for the exon scores is 0.5")
		#message("The used significance level for the p-values is 0.05")
		
		TestedData=ExonTesting(Data=Data,Exonthreshold=0,groups=groups,paired=paired,significancelevel=NULL)
		
		if(nrow(TestedData)!=0){
#			message(paste("Keep probesets with exon score greater than",Exonthreshold,sep=" "))
#			Data_Filt1=TestedData[which(TestedData$X50.>Exonthreshold),]
#		
#			message("Adjusting p-values for multiplicity")	
#			Data_Filt1$adj.p.value=p.adjust(Data_Filt1$p.value,"fdr")
#		
#			message(paste("Keep probesets with a p-value lower than",significancelevel,sep=" "))
#			Data_Sign=Data_Filt1[which(Data_Filt1$adj.p.value<0.05),]
			
			
			message("Ordering data from high to low significance")
			Data_Sign_Ordered=TestedData[order(-TestedData$adj.p.value),]
			rownames(Data_Sign_Ordered)=c(1:nrow(Data_Sign_Ordered))
		}
		else{
			Data_Sign_Ordered=NULL
		}
		
	}
	else if(class(Data)=="data.frame"){
		message("In using this function please make sure that the data has not been filtered yet and still has unadjusted p-values.")
		Data_Sign_Ordered=Data
		Data_Sign_Ordered$adj.p.value=p.adjust(Data_Sign_Ordered$p.value,"fdr")
		Data_Sign_Ordered=Data_Sign[order(-Data_Sign$X50.),]
		rownames(Data_Sign_Ordered)=c(1:nrow(Data_Sign_Ordered))
	}	
	Out=Data_Sign_Ordered
	
	
	message("Fusing Probe and Junction information")
#	
	colnames(JunctionInfoFile)[2]="ExonID"
	JoinedFile=merge(Out,JunctionInfoFile,by="ExonID",all=TRUE)
	JoinedFile=JoinedFile[,c(2,1,3,4,5,7,8,9,10,11,12)]
	colnames(JoinedFile)[1]="GeneID"
	
	return(JoinedFile)
	
}

##Testing of groups
groups=list(group1=groupLiver,group2=groupMuscle)


HTAData_RASA_PSR_LiverVSMuscle_AllPSR=AS_PSR_JUN(Data=HTAData_RASA_WithoutJunctions_OutputGDSModel,Exonthreshold=0,groups=groups,paired=FALSE,significancelevel=NULL,JunctionInfoFile=Probes_and_Junctions)
save(HTAData_RASA_PSR_LiverVSMuscle_AllPSR,file="HTAData_RASA_PSR_LiverVSMuscle_AllPSR.RData")









