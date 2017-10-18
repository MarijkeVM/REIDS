## dependent packages
#library(mvtnorm)
library(MCMCpack) # Iwish for udpating Sigma.b
library(nlme) #lme function
library(data.table)

#call on other functions
source("inigds.R")
source("gds.R")

#get data ready
args <- commandArgs(TRUE)
file_name <- as.character(args[1])
file_pos <- as.numeric(args[2])
line_length <- as.numeric(args[3])

conn <- file(file_name, 'rb')
current.pos <- seek(conn, where = file_pos, origin = 'start')
data <- readBin(conn, 'raw', n = line_length)
s <- paste(sapply(data, rawToChar), collapse='')

t=strsplit(s,"\",\"")

d=unlist(t)


d[1]=substr(d[1],start=2,stop=nchar(d[1]))
d[5]=substr(d[5],start=1,stop=nchar(d[5])-1)


if(d[1]!="geneID"){
	
	
	
	geneID=as.character(d[1])
	exonID=as.character(unlist(strsplit(d[2],",")))
	lengthexons=as.integer(unlist(strsplit(d[3],",")))
	
	npersample=sum(lengthexons)
	
	allsamples=d[4]
	samples=as.numeric(unlist(strsplit(allsamples,",")))
	samplenames=as.character(unlist(strsplit(d[5],",")))
	nsamples=length(samplenames)
	
	splitsamples<-function(x,samples,npersample){
		start=1+npersample*(x-1)
		end=npersample*x
		values=samples[start:end]
		return(values)
	}
	
	samplevalues=lapply(c(1:nsamples),function(i) splitsamples(i,samples,npersample) )
	TempData=rbindlist(list(samplevalues))
	setnames(TempData,colnames(TempData),samplenames)
	
	geneData=data.frame(geneID=rep(geneID,npersample),exonID=rep(exonID,lengthexons))
	geneData=data.frame(lapply(geneData, as.character), stringsAsFactors=FALSE)
	geneData=cbind(geneData,TempData)
	
	
	
	GenoDiffSplicing <- function(geneData,nsim=1000,geneID,informativeCalls=TRUE,alpha=0.5){
		
		Juncs=which(sapply(geneData[,2],function(x) substr(x,1,3))=="JUC")
		if(length(Juncs)>0){
			geneData=geneData[-c(which(sapply(geneData[,2],function(x) substr(x,1,3))=="JUC")),,drop=FALSE]
		}
		exonScore <- arrayScore <- informativeData<- NULL
		
		output=list()
		#for(e in 1:length(unique(geneID))){ # 
		output[[1]]=list()
		names(output)[1]=geneID
		lcmmData <- geneData[,-c(1,2)]  # relies on rownames: more secure in selecting correct rows
		lcmmData=as.matrix(lcmmData)
		enames <- geneData$exonID
		names(enames)=NULL
		
		rownames(lcmmData)<- enames
		
		##informative calls
		i=1
		if(informativeCalls){		

			fit <- inigds(SubgeneData=lcmmData, nsim) 
			fit2 <- data.frame(exonNames=unique(enames),Score=fit,informative=fit>alpha) # inigds returns one value per exon: filtering on exon level,no replicates
			output[[1]][[i]]=fit2
			names(output[[1]])[i]="Informative"
			#informativeData<- rbind(informativeData,fit2)
			iniData <- lcmmData[which(rownames(lcmmData)%in%fit2$exonNames[fit2$informative]),] # Of those that pass filtering step, retrieve the replicates and samples
			lcmmData <-  iniData   
			i=i+1
		}
		
		if(!is.null(lcmmData)&length(unique(rownames(lcmmData)))>1){  
			fit <- gds(SubgeneData=lcmmData, nsim) 
			#exonScore <- rbind(exonScore,fit$exonScores)
			exonScore <-fit$exonScores
			#arrayScore <- rbind(arrayScore ,fit$arrayScores)
			arrayScore <- fit$arrayScores
			output[[1]][[i]]=exonScore
			names(output[[1]])[i]="exonScore"
			output[[1]][[i+1]]=arrayScore
			names(output[[1]])[i+1]="arrayScore"
		}
		return(output)
	}
	
	Output=GenoDiffSplicing(geneData=geneData,nsim=5000,geneID=geneID,informativeCalls=TRUE,alpha=0.5)
	
	OutputREIDS_geneID=Output
	assign(paste("OutputREIDS_geneID",geneID,sep=""),OutputREIDS_geneID)
	eval(parse(text=paste("save(OutputREIDS_geneID", geneID, ", file=\"/scratch/leuven/310/vsc31011/HTA_RASA/REIDS/OutputREIDS/OutputREIDS_geneID", geneID, ".RData\")", sep=""))) 
	
	
}else{
	print(d)
}
