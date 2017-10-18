geneID=read.csv(file="geneID.csv",header=TRUE)

CreateOutput<-function(ID){
	Output=list()
	for(i in ID$GeneID){
		Data=get(load(paste("OutputGDS_geneID",as.character(i),".RData",sep="")))
		
		Output[length(Output)+1]=Data
		names(Output)[length(Output)]=i
	}
	HTAData_RASA_WithoutJunctions_OutputGDSModel=Output
	save(HTAData_RASA_WithoutJunctions_OutputGDSModel,file="../HTAData_RASA_WithoutJunctions_OutputREIDSModel.RData")
	
}

CreateOutput(geneID)





