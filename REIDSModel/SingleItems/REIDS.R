library(MCMCpack) # Iwish for udpating Sigma.b
library(nlme) #lme function
library(data.table)

#call on other functions
source("inigds.R")
source("gds.R")


REIDS <- function(geneData,nsim=1000,geneID,informativeCalls=TRUE,alpha=0.5){
	
	Juncs=which(sapply(geneData[,2],function(x) substr(x,1,3))=="JUC")
	if(length(Juncs)>0){
		geneData=geneData[-c(which(sapply(geneData[,2],function(x) substr(x,1,3))=="JUC")),,drop=FALSE]
	}
	exonScore <- arrayScore <- informativeData<- NULL
	
	output=list()
	for(e in 1:length(unique(geneID))){ 
		output[[e]]=list()
		names(output)[e]=unique(geneID)[e]
		lcmmData <- geneData[which(rownames(geneData)==unique(geneID)[e]),]  
		enames <- exonID[which(names(exonID)==unique(geneID)[e])]
		names(enames)=NULL
		rownames(lcmmData)<- enames
		
		i=1
		if(informativeCalls){			
			fit <- inigds(SubgeneData=lcmmData, nsim) 
			fit2 <- data.frame(exonNames=unique(enames),Score=fit,informative=fit>alpha) # iniREIDS returns one value per exon: filtering on exon level,no replicates
			output[[e]][[i]]=fit2
			names(output[[e]])[i]="Informative"
			iniData <- lcmmData[which(rownames(lcmmData)%in%fit2$exonNames[fit2$informative]),] # Of those that pass filtering step, retrieve the replicates and samples
			lcmmData <-  iniData   
			i=i+1
		}
	
		if(!is.null(lcmmData)&length(unique(rownames(lcmmData)))>1){  
			fit <- gds(SubgeneData=lcmmData, nsim) 
			exonScore <-fit$exonScores
			arrayScore <- fit$arrayScores
			output[[e]][[i]]=exonScore
			names(output[[e]])[i]="exonScore"
			output[[e]][[i+1]]=arrayScore
			names(output[[e]])[i+1]="arrayScore"
		}
	
	}
	
	return(output)
}

Output=GenoDiffSplicing(geneData=geneData,nsim=5000,geneID=geneID,informativeCalls=TRUE,alpha=0.5)

