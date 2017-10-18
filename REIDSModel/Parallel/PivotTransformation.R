
load("HTAData_RASA.rda")
library(data.table)

# HTA_RASA has exons names stiil with _s_st1, _s_st2,...,_s_st8
# Split these names to get regular exon IDs

GeneID=unique(as.character(HTAData_RASA$GeneID))

#xonID=as.character(HTAData_RASA$ExonID)
#xonID_temp=sapply(1:length(ExonID),function(x) strsplit(ExonID[x],"_")[[1]][1])
HTAData_RASA$ExonID=ExonID_temp
save(HTAData_RASA,file="HTAData_RASA.rda")
ExonID=unique(as.character(HTAData_RASA$ExonID))

PivotTransformData<-function(Data, geneID, exonID,savecsv=FALSE,name=NULL,location=NULL){
	Data=as.data.frame(Data)
#	if(is.null(Data$GeneID)&is.null(Data$ExonID)){
#		DataTemp=data.frame(geneID=geneID,exonID=exonID)
#		DataTemp=cbind(DataTemp,Data)
#		Data=DataTemp
#	}
	
	## From this point we assume that the first column of the data is a gene ID and the second column is the exonID
	
	Transformation<-function(x,gID,Data){
		
		if((x%%1000)==0){
			anotherthousand=c(x)
			write.table(anotherthousand,file=paste(x,".csv",sep=""))
		}
		Subset=Data[which(Data$GeneID==gID),]
		
		generow=c()
		gID=as.character(gID)
		
		eID=unique(Subset$ExonID)
		lengthe=sapply(eID,function(x) return(length(which(Subset$ExonID==x))))
		
		eID=paste(as.character(eID),sep="",collapse=",")
		lengthe=paste(as.character(lengthe),sep="",collapse=",")
		
		#get all samples at once: apply on columns of the Subset data and discard the geneID and exonID columns
		samples=apply(Subset[,-c(1,2)],2,function(x) paste(as.character(x),sep="",collapse=","))
		allsamples=paste(samples,sep="",collapse=",")
		samplenames=paste(colnames(Data)[-c(1,2)],sep="",collapse=",")
		#put everything intro c()
		generow=c(gID,eID,lengthe,allsamples,samplenames)
		
	}
	#use rbindlist to get a full data file
	DataPivot=lapply(1:length(unique(geneID)),function(x) Transformation(x,gID=unique(geneID)[x],Data=Data))
	DataBind=t(rbindlist(list(DataPivot)))
	colnames(DataBind)=c("geneID","exonID","lengthexons","allsamples","samplenames")
	rownames(DataBind)=unique(geneID)
	
	DataBind=as.data.frame(DataBind)
	DataBind$geneID=as.character(DataBind$geneID)
	DataBind$exonID=as.character(DataBind$exonID)
	DataBind$lengthexons=as.character(DataBind$lengthexons)
	DataBind$allsamples=as.character(DataBind$allsamples)
	DataBind$samplenames=as.character(DataBind$samplenames)
	
	if(savecsv){
		
		#location1=paste(location,".RData",sep="")
		location2=paste(location,".csv",sep="")
		
		assign(name,DataBind)
		#save(name,file=location1)
		write.table(get(name),file=location2,row.names=FALSE,col.names=TRUE,sep=",",quote=TRUE,qmethod="double")
	}
	
	return(DataBind)
}

HTAData_RASA_Cluster=PivotTransformData(Data=HTAData_RASA,geneID=GeneID,exonID=ExonID,savecsv=TRUE,name="HTAData_RASA_Cluster",location="HTAData_RASA_Cluster")

GeneTable=data.frame("GeneID"=GeneID)
head(GeneTable)
#      GeneID
# 1 TC0100001
# 2 TC0100002
# 3 TC0100003
# 4 TC0100006
# 5 TC0100007
# 6 TC0100009

write.table(GeneTable,file="geneID.csv",row.names=FALSE)