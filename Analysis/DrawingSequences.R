# Project: Detection of Splicing Variants
# 
# Author: lucp8409
###############################################################################


library(Biostrings)
library(plotrix)

Positions=read.table("ProbsetIDS_outcsv.txt",header=TRUE)
Mapping=read.table("F:/Doctorate/Detection of Splicing Variants/Preprocessing/Aroma Package/HTA 6 samples - run 2/Creating CDF/TC_PSR_JUCID_NumberID.txt",sep="\t",header=FALSE)
ASStructure=read.table("F:/Doctorate/Detection of Splicing Variants/Preprocessing/Aroma Package/HTA 6 samples - run 2/Creating CDF/hta1.r1.ass",sep="\t",header=TRUE)
ChrStructure=read.table("F:/Doctorate/Detection of Splicing Variants/Preprocessing/Aroma Package/HTA 6 samples - run 2/Creating CDF/hta1.r1.psr_ps.bed",sep="\t",header=FALSE)
colnames(ChrStructure)=c("chr","start","end","probesetid","offset","strand")


DrawingSequences<-function(Junction,ASStucture,MappingFile=Mapping,SequenceFile=Positions,ChrFile=ChrStructure){
	JunctionMatchingSequences=c()
	
	Exons=ASStructure[which(ASStructure$junction_id==Junction&ASStructure$as_type!="exclusion"),]
	Exon3=Exons[which(Exons$as_type==3),2]
	Exon5=Exons[which(Exons$as_type==5),2]
	if(length(Exon3)==0){
		Exon3=NULL
	}
	if(length(Exon5)==0){
		Exon5=NULL
	}
	
	if(length(Exon3)!=0){
		Exon3_ID=MappingFile[which(as.character(MappingFile[,2])==as.character(Exon3)),1]
	}	
	else{
		Exon3_ID="not present"
	}
	Junction_ID=MappingFile[which(as.character(MappingFile[,2])==as.character(Junction)),1]
	if(length(Exon5)!=0){
		Exon5_ID=MappingFile[which(as.character(MappingFile[,2])==as.character(Exon5)),1]
	}
	else{
		Exon5_ID="not present"
	}
	
	Sequences=list()
	IDs=list(as.character(Exon3),Junction,as.character(Exon5))
	for(e in 1:length(IDs)){
		print(e)
		if(length(IDs[[e]])>0){
			Mapped=SequenceFile[which(as.character(SequenceFile[,1])==IDs[e]),]
			t1=read.table(file="out.csv",sep="\t",header=FALSE,colClasses="character",skip=as.numeric((Mapped[2]-1)),nrows=as.numeric(Mapped[3]))
			Sequences[[e]]=t1[,c(1,5,6,11)]
		}
			
	}
	names(Sequences)=names(IDs)
	
	Aligned=list()
	for(k in 1:length(Sequences)){
		if(!is.null(Sequences[[k]])){
			Temp=Sequences[[k]]
			Temp$exon_position=as.numeric(Temp[,2])
			Temp=Temp[order(Temp$exon_position),]
			for(i in 1:nrow(Temp)){
				if(i==1){
					Seq=Temp[i,4]
				}
				else{
					ExonPosEndPrev=as.numeric(Temp[i-1,2])+nchar(Temp[i-1,4])
					ExonPosEnd=as.numeric(Temp[i,2])+nchar(Temp[i,4])
					Rest=ExonPosEnd-ExonPosEndPrev
					if(Rest<25){
						Seq=paste(Seq,substr(Temp[i,4],start=(nchar(Temp[i,4])-Rest+1),stop=nchar(Temp[i,4])),sep="")
					}
					else{
						Seq=paste(Seq,Temp[i,4],sep="...")
					}
					
				}
			}
			Aligned[[k]]=Seq
		}		
		else{
			Aligned[[k]]=0
		}
	}
	names(Aligned)=names(Sequences)
	
	#strand=ChrFile[which(ChrFile[,4]==Exon5_ID),6]
	out=c()
	if(Aligned[[1]]!=0){
		FirstPartJunction=substr(Aligned[[1]],start=nchar(Aligned[[1]])-16,stop=nchar(Aligned[[1]]))
		Aligned1=pairwiseAlignment(FirstPartJunction,substr(Aligned[[2]],1,nchar(Aligned[[2]])/2),type="local",gapOpening=10, gapExtension=4)@pattern
		if(all(strsplit(substr(FirstPartJunction,start=nchar(FirstPartJunction)-nchar(Aligned1)+1,stop=nchar(FirstPartJunction)),'')[[1]]==strsplit(as.character(Aligned1),'')[[1]])){
			#message("Junctions matches the first exon")
			out=c(out,as.character(Exon3),"match")
		}
		else{
			out=c(out,as.character(Exon3),"no match")
		}
	}
	else{
		out=c(out,"-","-")
	}
	
	
	if(Aligned[[3]]!=0){
		LastPartJunction=substr(Aligned[[3]],start=1,stop=25)
		Aligned2=pairwiseAlignment(LastPartJunction,substr(Aligned[[2]],nchar(Aligned[[2]])/2,nchar(Aligned[[2]])),type="local",gapOpening=10, gapExtension=4)@pattern
		if(all(strsplit(substr(LastPartJunction,start=1,stop=nchar(Aligned2)),'')[[1]]==strsplit(as.character(Aligned2),'')[[1]])){
			#message("Junctions matches the second exon")
			out=c(out,as.character(Exon5),"match")
		}
		else{
			out=c(out,as.character(Exon5),"no match")
		}
	}
	else{
		out=c(out,"-","-")
	}
	
	JunctionMatchingSequences=rbind(JunctionMatchingSequences,c(as.character(Junction),out))
	colnames(JunctionMatchingSequences)=c("Junction","Exon3","State","Exon5","State")
	JunctionMatchingSequences=data.frame(JunctionMatchingSequences)
	NR=sum(unlist(lapply(Sequences,nrow)))
	NC=sum(unlist(lapply(Aligned,nchar)))
	
	Overlap1=Overlap2=0
	EndOverlap1=EndOverlap2=0
	M=matrix(0,nrow=NR,ncol=NC)
	
	Null=which(sapply(Sequences,function(x) is.null(x)))
	if(length(Null)>0){
		Sequences=Sequences[-c(Null)]
	}
	
	Jump1=0
	Jump2=0
	
	for(i in 1:length(Sequences)){
		print(i)
		Temp=Sequences[[i]]
		Temp$exon_position=as.numeric(Temp[,2])
		Temp=Temp[order(Temp$exon_position),]
		Temp_probes=sapply(Temp[,4],strsplit, "")

		for(j in 1:length(Temp_probes)){
			if(i==1){
				if(j==1){
					M[1,c(1:25)]=Temp_probes[[1]]
					Begin=1
					End=25
				}
				else{
					Start=Temp$exon_position[j-1]
					Shift=Temp$exon_position[j]
					if((Shift-Start)>25){
						M[j,c((Begin+1):(Begin+26))]=c("...",Temp_probes[[j]])
						if(j>Jump1){
							Jump1=j
						}
						Begin=Begin+2
						End=Begin+24
					}
					else{
						Move=Shift-Start
						M[j,c((Begin+Move):(Begin+Move+24))]=Temp_probes[[j]]
						Begin=Begin+Move
						End=Begin+24
					}
				}
			}
			else if(i==2 & JunctionMatchingSequences[which(as.character(JunctionMatchingSequences[,1])==as.character(Junction)),3]!="-"){
				if(j==1 & JunctionMatchingSequences[which(JunctionMatchingSequences[,1]==Junction),3]=="match"){
					GoBack=nchar(Aligned1)
					Begin=End-GoBack
					Overlap1=Begin+1
					EndOverlap1=Overlap1+nchar(Aligned1)-1
				}	
				else if(j==1){
					Begin=End
					print(Begin)
				}	
				if(j==1){
					Begin=Begin+1
					M[j+nrow(Sequences[[1]]),c((Begin):(Begin+24))]=Temp_probes[[j]]
					End=Begin+24
				}
				else{
					Begin=Begin+1
					M[j+nrow(Sequences[[1]]),c((Begin):(Begin+24))]=Temp_probes[[j]]
					End=Begin+24				
				}
			}
			else if(i==2 & JunctionMatchingSequences[which(JunctionMatchingSequences[,1]==Junction),3]=="-"){
				if(j==1 & JunctionMatchingSequences[which(JunctionMatchingSequences[,1]==Junction),5]=="match"){
					GoBack=nchar(Aligned2)
					Begin=End-GoBack
					Overlap2=Begin+1
					EndOverlap2=Overlap2+nchar(Aligned2)-1
				}
				else if(j==1){
					Begin=End
				}
				if(j==1){
					Begin=Begin+1
					
					M[j+nrow(Sequences[[1]]),c(Begin:c(Begin+24))]=Temp_probes[[j]]
					End=Begin+24
				}
				else{
					Start=Temp$exon_position[j-1]
					Shift=Temp$exon_position[j]
					if((Shift-Start)>25){
						M[j+nrow(Sequences[[1]]),c((Begin+1):(Begin+26))]=c("...",Temp_probes[[j]])
						Begin=Begin+2
						End=Begin+24
					}
					else{
						Move=Shift-Start
						M[j+nrow(Sequences[[1]]),c((Begin+Move):(Begin+Move+24))]=Temp_probes[[j]]
						Begin=Begin+Move
						End=Begin+24
					}
			}
			}
			else if(i==3){
				if(j==1 & JunctionMatchingSequences[which(JunctionMatchingSequences$Junction==Junction),5]=="match"){
					GoBack=nchar(Aligned2)
					Begin=End-GoBack
					Overlap2=Begin+1
					EndOverlap2=Overlap2+nchar(Aligned2)-1
				}
				else if(j==1){
					Begin=End
				}
				if(j==1){
					Begin=Begin+1
					M[j+nrow(Sequences[[1]])+nrow(Sequences[[2]]),c(Begin:c(Begin+24))]=Temp_probes[[j]]
					End=Begin+24
				}
				else{
					Start=Temp$exon_position[j-1]
					Shift=Temp$exon_position[j]
					if((Shift-Start)>25){
						M[j+nrow(Sequences[[1]])+nrow(Sequences[[2]]),c((Begin+1):(Begin+26))]=c("...",Temp_probes[[j]])
						if(j+nrow(Sequences[[1]])+nrow(Sequences[[2]])>Jump2){
							Jump2=j+nrow(Sequences[[1]])+nrow(Sequences[[2]])
						}
						Begin=Begin+2
						End=Begin+24
					}
					else{
						Move=Shift-Start
						M[j+nrow(Sequences[[1]])+nrow(Sequences[[2]]),c((Begin+Move):(Begin+Move+24))]=Temp_probes[[j]]
						Begin=Begin+Move
						End=Begin+24
					}
				}
			}
		}		
	}
	M2<-M
	M2[]<-match(M,c("A","C","G","T",0))
	mode(M2)<-"numeric"
	M2[which(M2==5)]=0
	DelCols=which(colSums(M2)==0)
	DelRows=which(rowSums(M2)==0)
	M2=M2[,-DelCols]
#	cellcolors=as.vector(M2)
#	cellcolors[which(cellcolors==0)]="white"
#	color2D.matplot(M2,cellcolors=cellcolors,border=NA)
	
	color2D.matplot(M2,cellcolors=rep("white",length(as.vector(M))),border=NA)
	
	M=M[,-DelCols]
	M[which(M=="0")]=""
	for(i in 1:nrow(M)){
		if(i<=nrow(Sequences[[1]])){
			col=rep("red",ncol(M))
			if(Overlap1!=0 & i>=Jump1){
				col[c(Overlap1:EndOverlap1)]="purple"
			}	
		}
		else if(i>(nrow(Sequences[[1]])+nrow(Sequences[[2]]))){
			col=rep("dark green",ncol(M))	
			if(Overlap2!=0 & i<=Jump2){
				col[c(Overlap2:EndOverlap2)]="orange"
			}
			
		}
		else{
			col=rep("blue",ncol(M))
			if(Overlap2!=0){
				col[c(Overlap2:EndOverlap2)]="orange"
			}
			
			if(Overlap1!=0){
				col[c(Overlap1:EndOverlap1)]="purple"
			}
			
		}
		text(c(0.5:(ncol(M)-0.5)),rep(c((nrow(M)+0.5)-i),ncol(M)),M[i,],col=col,cex=1)
	}
}

DrawingSequences_5no3<-function(Junction,ASStucture,MappingFile=Mapping,SequenceFile=Positions,ChrFile=ChrStructure){
	JunctionMatchingSequences=c()
	
	Exons=ASStructure[which(ASStructure$junction_id==Junction&ASStructure$as_type!="exclusion"),]
	Exon3=Exons[which(Exons$as_type==3),2]
	Exon5=Exons[which(Exons$as_type==5),2]
	if(length(Exon3)==0){
		Exon3=NULL
	}
	if(length(Exon5)==0){
		Exon5=NULL
	}
	
	if(length(Exon3)!=0){
		Exon3_ID=MappingFile[which(as.character(MappingFile[,2])==as.character(Exon3)),1]
	}	
	else{
		Exon3_ID="not present"
	}
	Junction_ID=MappingFile[which(as.character(MappingFile[,2])==as.character(Junction)),1]
	if(length(Exon5)!=0){
		Exon5_ID=MappingFile[which(as.character(MappingFile[,2])==as.character(Exon5)),1]
	}
	else{
		Exon5_ID="not present"
	}
	
	Sequences=list()
	IDs=list(as.character(Exon3),Junction,as.character(Exon5))
	for(e in 1:length(IDs)){
		print(e)
		if(length(IDs[[e]])>0){
			Mapped=SequenceFile[which(as.character(SequenceFile[,1])==IDs[e]),]
			t1=read.table(file="out.csv",sep="\t",header=FALSE,colClasses="character",skip=as.numeric((Mapped[2]-1)),nrows=as.numeric(Mapped[3]))
			Sequences[[e]]=t1[,c(1,5,6,11)]
		}
		
	}
	names(Sequences)=names(IDs)
	
	Aligned=list()
	for(k in 1:length(Sequences)){
		if(!is.null(Sequences[[k]])){
			Temp=Sequences[[k]]
			Temp$exon_position=as.numeric(Temp[,2])
			Temp=Temp[order(Temp$exon_position),]
			for(i in 1:nrow(Temp)){
				if(i==1){
					Seq=Temp[i,4]
				}
				else{
					ExonPosEndPrev=as.numeric(Temp[i-1,2])+nchar(Temp[i-1,4])
					ExonPosEnd=as.numeric(Temp[i,2])+nchar(Temp[i,4])
					Rest=ExonPosEnd-ExonPosEndPrev
					if(Rest<25){
						Seq=paste(Seq,substr(Temp[i,4],start=(nchar(Temp[i,4])-Rest+1),stop=nchar(Temp[i,4])),sep="")
					}
					else{
						Seq=paste(Seq,Temp[i,4],sep="...")
					}
					
				}
			}
			Aligned[[k]]=Seq
		}		
		else{
			Aligned[[k]]=0
		}
	}
	names(Aligned)=names(Sequences)
	
	#strand=ChrFile[which(ChrFile[,4]==Exon5_ID),6]
	out=c()
	if(Aligned[[1]]!=0){
		FirstPartJunction=substr(Aligned[[1]],start=nchar(Aligned[[1]])-16,stop=nchar(Aligned[[1]]))
		Aligned1=pairwiseAlignment(FirstPartJunction,substr(Aligned[[2]],1,nchar(Aligned[[2]])/2),type="local",gapOpening=10, gapExtension=4)@pattern
		if(all(strsplit(substr(FirstPartJunction,start=nchar(FirstPartJunction)-nchar(Aligned1)+1,stop=nchar(FirstPartJunction)),'')[[1]]==strsplit(as.character(Aligned1),'')[[1]])){
			#message("Junctions matches the first exon")
			out=c(out,as.character(Exon3),"match")
		}
		else{
			out=c(out,as.character(Exon3),"no match")
		}
	}
	else{
		out=c(out,"-","-")
	}
	
	
	if(Aligned[[3]]!=0){
		LastPartJunction=substr(Aligned[[3]],start=1,stop=25)
		Aligned2=pairwiseAlignment(LastPartJunction,Aligned[[2]],type="local",gapOpening=10, gapExtension=4)@pattern
		if(all(strsplit(substr(LastPartJunction,start=1,stop=nchar(Aligned2)),'')[[1]]==strsplit(as.character(Aligned2),'')[[1]])){
			#message("Junctions matches the second exon")
			out=c(out,as.character(Exon5),"match")
		}
		else{
			out=c(out,as.character(Exon5),"no match")
		}
	}
	else{
		out=c(out,"-","-")
	}
	
	JunctionMatchingSequences=rbind(JunctionMatchingSequences,c(as.character(Junction),out))
	colnames(JunctionMatchingSequences)=c("Junction","Exon3","State","Exon5","State")
	JunctionMatchingSequences=data.frame(JunctionMatchingSequences)
	NR=sum(unlist(lapply(Sequences,nrow)))
	NC=sum(unlist(lapply(Aligned,nchar)))
	
	Overlap1=Overlap2=0
	EndOverlap1=EndOverlap2=0
	M=matrix(0,nrow=NR,ncol=NC)
	
	Null=which(sapply(Sequences,function(x) is.null(x)))
	if(length(Null)>0){
		Sequences=Sequences[-c(Null)]
	}
	
	for(i in 1:length(Sequences)){
		print(i)
		Temp=Sequences[[i]]
		Temp$exon_position=as.numeric(Temp[,2])
		Temp=Temp[order(Temp$exon_position),]
		Temp_probes=sapply(Temp[,4],strsplit, "")
		
		for(j in 1:length(Temp_probes)){
			if(i==1){
				if(j==1){
					M[1,c(1:25)]=Temp_probes[[1]]
					Begin=1
					End=25
				}
				else{
					Start=Temp$exon_position[j-1]
					Shift=Temp$exon_position[j]
					if((Shift-Start)>25){
						M[j,c((Begin+1):(Begin+26))]=c("...",Temp_probes[[j]])
						Begin=Begin+2
						End=Begin+24
					}
					else{
						Move=Shift-Start
						M[j,c((Begin+Move):(Begin+Move+24))]=Temp_probes[[j]]
						Begin=Begin+Move
						End=Begin+24
					}
				}
			}
			else if(i==2 & JunctionMatchingSequences[which(as.character(JunctionMatchingSequences[,1])==as.character(Junction)),3]!="-"){
				if(j==1 & JunctionMatchingSequences[which(JunctionMatchingSequences[,1]==Junction),3]=="match"){
					GoBack=nchar(Aligned1)
					Begin=End-GoBack
					Overlap1=Begin+1
					EndOverlap1=Overlap1+nchar(Aligned1)-1
				}	
				else if(j==1){
					Begin=End
					print(Begin)
				}	
				if(j==1){
					Begin=Begin+1
					M[j+nrow(Sequences[[1]]),c((Begin):(Begin+24))]=Temp_probes[[j]]
					End=Begin+24
				}
				else{
					Begin=Begin+1
					M[j+nrow(Sequences[[1]]),c((Begin):(Begin+24))]=Temp_probes[[j]]
					End=Begin+24				
				}
			}
			else if(i==2 & JunctionMatchingSequences[which(JunctionMatchingSequences[,1]==Junction),3]=="-"){
				if(j==1 & JunctionMatchingSequences[which(JunctionMatchingSequences[,1]==Junction),5]=="match"){
					GoBack=nchar(Aligned2)
					Begin=End-GoBack
					Overlap2=Begin+1
					EndOverlap2=Overlap2+nchar(Aligned2)-1
				}
				else if(j==1){
					Begin=End
				}
				if(j==1){
					Begin=Begin+1
					
					M[j+nrow(Sequences[[1]]),c(Begin:c(Begin+24))]=Temp_probes[[j]]
					End=Begin+24
				}
				else{
					Start=Temp$exon_position[j-1]
					Shift=Temp$exon_position[j]
					if((Shift-Start)>25){
						M[j+nrow(Sequences[[1]]),c((Begin+1):(Begin+26))]=c("...",Temp_probes[[j]])
						Begin=Begin+2
						End=Begin+24
					}
					else{
						Move=Shift-Start
						M[j+nrow(Sequences[[1]]),c((Begin+Move):(Begin+Move+24))]=Temp_probes[[j]]
						Begin=Begin+Move
						End=Begin+24
					}
				}
			}
			else if(i==3){
				if(j==1 & JunctionMatchingSequences[which(JunctionMatchingSequences$Junction==Junction),5]=="match"){
					GoBack=nchar(Aligned2)
					Begin=End-GoBack
					Overlap2=Begin+1
					EndOverlap2=Overlap2+nchar(Aligned2)-1
				}
				else if(j==1){
					Begin=End
				}
				if(j==1){
					Begin=Begin+1
					M[j+nrow(Sequences[[1]])+nrow(Sequences[[2]]),c(Begin:c(Begin+24))]=Temp_probes[[j]]
					End=Begin+24
				}
				else{
					Start=Temp$exon_position[j-1]
					Shift=Temp$exon_position[j]
					if((Shift-Start)>25){
						M[j+nrow(Sequences[[1]])+nrow(Sequences[[2]]),c((Begin+1):(Begin+26))]=c("...",Temp_probes[[j]])
						Begin=Begin+2
						End=Begin+24
					}
					else{
						Move=Shift-Start
						M[j+nrow(Sequences[[1]])+nrow(Sequences[[2]]),c((Begin+Move):(Begin+Move+24))]=Temp_probes[[j]]
						Begin=Begin+Move
						End=Begin+24
					}
				}
			}
		}		
	}
	M2<-M
	M2[]<-match(M,c("A","C","G","T",0))
	mode(M2)<-"numeric"
	M2[which(is.na(M2))]=0
	DelCols=which(colSums(M2)==0)
	DelRows=which(rowSums(M2)==0)
	if(length(DelCols)>0){
		M2=M2[,-DelCols]
	}
#	cellcolors=as.vector(M2)
#	cellcolors[which(cellcolors==0)]="white"
#	color2D.matplot(M2,cellcolors=cellcolors,border=NA)
	
	color2D.matplot(M2,cellcolors=rep("white",length(as.vector(M))),border=NA)
	
	if(length(DelCols)>0){
		M=M[,-DelCols]
	}
	M[which(M=="0")]=""
	for(i in 1:nrow(M)){
		if(i<=nrow(Sequences[[1]])){
			col=rep("blue",ncol(M))
			if(Overlap2!=0){
				col[c(Overlap2:EndOverlap2)]="orange"
			}	
		}
		else if(i>(nrow(Sequences[[1]])+nrow(Sequences[[2]]))){
			col=rep("dark green",ncol(M))	
			if(Overlap2!=0){
				col[c(Overlap2:EndOverlap2)]="orange"
			}
			
		}
		else{
			col=rep("dark green",ncol(M))
			if(Overlap2!=0){
				col[c(Overlap2:EndOverlap2)]="orange"
			}
			
			if(Overlap1!=0){
				col[c(Overlap1:EndOverlap1)]="purple"
			}
			
		}
		text(c(0.5:(ncol(M)-0.5)),rep(c((nrow(M)+0.5)-i),ncol(M)),M[i,],col=col,cex=1)
	}
}


pdf("Figures/Sequence_JUC0100236357.pdf",width=21,height=7)
DrawingSequences(Junction="JUC0100236357")
dev.off()

pdf("Figures/Sequence_JUC0100236628.pdf",width=21,height=7)
DrawingSequences(Junction="JUC0100236628")
dev.off()

pdf("Figures/Sequence_JUC0100236573.pdf",width=21,height=7)
DrawingSequences(Junction="JUC0100236573")
dev.off()


pdf("Figures/Sequence_JUC0100236756.pdf",width=21,height=7)
DrawingSequences(Junction="JUC0100236756")
dev.off()

pdf("Figures/Sequence_JUC0100236756.pdf",width=21,height=7)
DrawingSequences(Junction="JUC0800063802")
dev.off()

#
pdf("Figures/Sequence_JUC0800064097.pdf",width=21,height=7)
DrawingSequences_5no3(Junction="JUC0800064097")
dev.off()


pdf("Figures/Sequence_JUC1900096711.pdf",width=21,height=7)
DrawingSequences(Junction="JUC1900096711")
dev.off()

pdf("Figures/Sequence_JUC1900096693.pdf",width=14,height=7)
DrawingSequences_5no3(Junction="JUC1900096693")
dev.off()

pdf("Figures/Sequence_JUC1200112702.pdf",width=21,height=7)
DrawingSequences(Junction="JUC1200112702")
dev.off()

pdf("Figures/Sequence_JUC1900096711_JUC1900096693.pdf",width=21,height=14)
par(mfrow=c(2,1))
DrawingSequences(Junction="JUC1900096711")
DrawingSequences_5no3(Junction="JUC1900096693")
dev.off()


pdf("Figures/Sequence_JUC1200112431.pdf",width=21,height=7)
DrawingSequences(Junction="JUC1200112431")
dev.off()