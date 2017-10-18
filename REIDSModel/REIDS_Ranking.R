
library(lmtest)

#Data
load("Data/HTAData_RASA_PSR_LiverVSMuscle_AllPSR.RData")
load("Data/HTAData_RASA.rda")
for(i in 3:ncol(HTAData_RASA)){
	HTAData_RASA[,i]=as.numeric(as.character(HTAData_RASA[,i]))
}


# Number of AS Probe sets:

dim(HTAData_RASA_PSR_LiverVSMuscle_AllPSR)
#filter PSR on ICC Threshold:
Filter=HTAData_RASA_PSR_LiverVSMuscle_AllPSR[HTAData_RASA_PSR_LiverVSMuscle_AllPSR$X50.>=0.50,]
#filter PSR  on significance level
ASPSR_PSR=Filter[Filter$adj.p.value<0.05,]
length(unique(ASPSR_PSR$ExonID))
# [1] 20992


options(warn=2)
ScoringBasedOnRanking<-function(ASProbeSets=unique(ASPSR_PSRwithJUN$PSR_ID), AnnotData=HTAData_RASA_PSR_LiverVSMuscle_AllPSR,Data=HTAData_RASA){
	
	Scores<-function(x,ASPSR,Annot,AnnotData,Data){
		print(x)
		if(any(Annot$as_type%in%c("3","5"))){
			if(all(c("3","5")%in%Annot$as_type)){ #both junctions are present
				if(any(Annot$as_type=="exclusion")){
					Annot=Annot[-which(Annot$as_type=="exclusion"),]
				}
				if(length(which(Annot$as_type%in%c("3","5")))>2){ #more than 1 to both or either sides: take clostest one to each side
					
					JUCLeft=NULL
					JUCRight=NULL
					if(length(which(Annot$as_type=="3"))==1){
						JUCRight=Annot[which(Annot$as_type=="3"),1]
						Annot2=Annot[-which(Annot[,1]==JUCRight),]
						JUCLeft=NULL
					}
					else if(length(which(Annot$as_type=="5"))==1){
						JUCLeft=Annot[which(Annot$as_type=="5"),1]
						Annot2=Annot[-which(Annot[,1]==JUCLeft),]
						JUCRight=NULL
					}
					else{
						Annot2=Annot
					}
					
					OtherAnnot=c()
					for(j in Annot2$JUC_ID){
						OtherAnnot=rbind(OtherAnnot,AnnotData[which(AnnotData$JUC_ID==j),c(2,8,15)])
					}
					OtherAnnot=OtherAnnot[-which(OtherAnnot$PSR_ID==as.character(ASPSR)),]
					if(any(OtherAnnot$as_type=="exclusion")){
						ExcluTable=table(OtherAnnot[,2])
						Mult=names(ExcluTable)[which(ExcluTable>1)]
						for(e in Mult){
							OtherAnnot=OtherAnnot[-which(OtherAnnot[,2]==e & OtherAnnot$as_type=="exclusion"),]
						}
					}
					if(any(table(OtherAnnot[,1])>1)){
						Mult=names(table(OtherAnnot[,1]))[which(table(OtherAnnot[,1])>1)]
						for(e in Mult){
							while(nrow(OtherAnnot[which(OtherAnnot[,1]==e),])>1){
								Del=which(OtherAnnot[,1]==e & OtherAnnot$as_type=="exclusion")
								OtherAnnot=OtherAnnot[-Del[1],]
							}
						}
					}
					
#					if(nrow(OtherAnnot)==0){
#						Row=c(as.character(ASPSR),"","",100000,"No second annotation for both junctions")
#						return(Row)
#					}
#					else if(nrow(OtherAnnot)==1){
#						Row=c(as.character(ASPSR),"","",100000,"No second annotation for one junction")
#						return(Row)
#					}			
					#else{
					ConvertNum=as.integer(substr(OtherAnnot[,1],6,nchar(OtherAnnot[,1])))-as.integer(substr(as.character(ASPSR),6,nchar(as.character(ASPSR))))
					names(ConvertNum)=OtherAnnot[,1]
					
					Left=ConvertNum[which(ConvertNum<0)]
					ClosestLeft=names(Left)[which.min(abs(Left))]
					
					Right=ConvertNum[which(ConvertNum>0)]
					ClosestRight=names(Right)[which.min(Right)]
					
					
					OtherAnnot=OtherAnnot[order(as.numeric(OtherAnnot$as_type)),]
					Junctions=OtherAnnot[OtherAnnot$PSR_ID%in%c(ClosestLeft,ClosestRight),2]
					if(length(Junctions==2)){
						if(length(unique(Annot[which(Annot[,1]%in%Junctions),2]))==1){
							Junctions=Junctions[1]
							if(any(!is.null(c(JUCRight,JUCLeft)))){
								JUCS=c(JUCRight,JUCLeft)
								Junctions=c(Junctions,JUCS)
							}
						}
					}
					if(length(ClosestRight)==0){
						if(is.null(JUCLeft)& !(is.null(JUCRight))){
							Junctions=c(Junctions,JUCRight)
						}
						else if (!(is.null(JUCLeft)) & is.null(JUCRight)){
							Junctions=c(Junctions,JUCLeft)
						}
					}
					else if(length(ClosestLeft)==0){
						if(is.null(JUCLeft) & !(is.null(JUCRight))){					
							Junctions=c(Junctions,JUCRight)
						}
						
						else if (!(is.null(JUCLeft)) & is.null(JUCRight)){
							Junctions=c(Junctions,JUCLeft)
						}
					}
					#print(Junctions)
					
					if(length(Junctions)==0){
						Juc5=Annot2[which(Annot2$as_type==5),1][1]
						Juc3=Annot2[which(Annot2$as_type==3),1][1]
						Junctions=c(Juc5,Juc3)
					}
					else if(length(Junctions)==1){
						if(Annot[which(Annot[,1]==Junctions),2]==3){
							Juc5=Annot2[which(Annot2$as_type==5),1][1]
							Junctions=c(Junctions,Juc5)
						}
						else if(Annot[which(Annot[,1]==Junctions),2]==5){
							Juc3=Annot2[which(Annot2$as_type==3),1][1]
							Junctions=c(Junctions,Juc3)
						}
					}
					
					Del=c()
					for(j in Junctions){
						State=AnnotData[which(AnnotData$JUC_ID==j),ncol(AnnotData)][1]
#						if(State=="filtered out by INI"){
#							Del=c(Del,j)
#						}
					}
					if(length(Del)>0){
						if(length(Del)==length(Junctions)){
							Row=c(as.character(ASPSR),Junctions[1],Junctions[2],100000,"INI Filtered Junctions")
							return(Row)
						}
						else{
							Row=c(as.character(ASPSR),Junctions[1],Junctions[2],100000,"INI Filtered Junction")
							return(Row)
							
						}
					}
					
					
					#}
					
				}
				else{
					Annot=Annot[order(as.numeric(Annot$as_type)),]
					Junctions=Annot$JUC_ID[which(Annot$as_type%in%c("3","5"))]
					Del=c()
					for(j in Junctions){
						State=AnnotData[which(AnnotData$JUC_ID==j),ncol(AnnotData)][1]
#						if(State=="filtered out by INI"){
#							Del=c(Del,j)
#						}
					}
					if(length(Del)>0){
						if(length(Del)==length(Junctions)){
							Row=c(as.character(ASPSR),Junctions[1],Junctions[2],100000,"INI Filtered Junctions")
							return(Row)
						}
						else{
							Junctions=Junctions[-which(Junctions%in%Del)]
							Row=c(as.character(ASPSR),Junctions[1],Junctions[2],100000,"INI Filtered Junction")
							return(Row)
						}
					}
				}
				
				Juc3=Annot$JUC_ID[which(Annot$as_type=="3"&Annot$JUC_ID%in%Junctions)]
				Juc5=Annot$JUC_ID[which(Annot$as_type=="5"&Annot$JUC_ID%in%Junctions)]	
				
				PSR=as.vector(as.matrix(Data[which(Data$ExonID==as.character(ASPSR)),-c(1,2)]))
				PSR_Ranks=sort(PSR,index.return=TRUE)$ix
				JUC3=as.vector(as.matrix(Data[which(Data$ExonID==Juc3),-c(1,2)]))
				JUC3_Ranks=sort(JUC3,index.return=TRUE)$ix
				JUC5=as.vector(as.matrix(Data[which(Data$ExonID==Juc5),-c(1,2)]))
				JUC5_Ranks=sort(JUC5,index.return=TRUE)$ix
				
				if(!Juc3%in%unique(Data$ExonID)){
					Row=c(as.character(ASPSR),Juc3,Juc5,100000,100000,"Junction 3' not found in Data")
					return(Row)
				}
				
				if(!Juc5%in%unique(Data$ExonID)){
					Row=c(as.character(ASPSR),Juc3,Juc5,100000,100000,"Junction 5' not found in Data")
					return(Row)
				}
				
				L=list(PSR_Ranks,JUC3_Ranks,JUC5_Ranks)
				if(length(unique(c(length(PSR),length(JUC3),length(JUC5))))>1){
					MinLength=min(c(length(PSR),length(JUC3),length(JUC5)))
					Index=which(c(length(PSR),length(JUC3),length(JUC5))>MinLength)
					for(s in Index){
						L[[s]]=L[[s]][-which(L[[s]]%in%c((MinLength+1):length(L[[s]])))]
					}
				}
				
				Ranks=cbind(L[[1]],L[[2]],L[[3]])
				
				Y=as.vector(as.matrix(Ranks))
				exon=c(rep(1,nrow(Ranks)),rep(2,nrow(Ranks)),rep(3,nrow(Ranks)))
				tissue=rep(c(rep(1,nrow(Ranks)/2),rep(2,nrow(Ranks)/2)),3)
				
				exon=as.factor(exon)
				tissue=as.factor(tissue)
				ft1<-lm(Y~exon*tissue-1)
				ft2<-lm(Y~exon+tissue)
				
				Chisq=lrtest(ft2,ft1)[2,4]
				
				InteractionSign=anova(ft1)[3,5]
				
				if(InteractionSign<0.05){
					if(any(summary(ft1)$coefficients[c(5,6),4]>0.05)){
						Row=c(as.character(ASPSR),Juc3,Juc5,InteractionSign,Chisq,"1 of two junctions is supporting")
						return(Row)
					}
					else{
						Row=c(as.character(ASPSR),Juc3,Juc5,InteractionSign,Chisq,"No junction is supporting")
						return(Row)
					}
				}
				
				Row=c(as.character(ASPSR),Juc3,Juc5,InteractionSign,Chisq,"")
			}	
			else{
				if(any(Annot$as_type=="exclusion")){
					Annot=Annot[-which(Annot$as_type=="exclusion"),]
				}
				if(length(which(Annot$as_type%in%c("3","5")))>=2){ #more than 1 to both or either sides: take clostest one to each side
					
#					JUCLeft=NULL
#					JUCRight=NULL
#					if(length(which(Annot$as_type=="3"))<1){
#						JUCRight=Annot[which(Annot$as_type=="5"),1]
#						Annot2=Annot[-which(Annot[,1]==JUCRight),]
#						JUCLeft=NULL
#					}
#					else if(length(which(Annot$as_type=="5"))<1){
#						JUCLeft=Annot[which(Annot$as_type=="3"),1]
#						Annot2=Annot[-which(Annot[,1]==JUCLeft),]
#						JUCRight=NULL
#					}
#					else{
					Annot2=Annot
#					}
					
					OtherAnnot=c()
					for(j in Annot2$JUC_ID){
						OtherAnnot=rbind(OtherAnnot,AnnotData[which(AnnotData$JUC_ID==j),c(2,8,15)])
					}
					OtherAnnot=OtherAnnot[-which(OtherAnnot$PSR_ID==as.character(ASPSR)),]
					if(any(OtherAnnot$as_type=="exclusion")){
						ExcluTable=table(OtherAnnot[,2])
						Mult=names(ExcluTable)[which(ExcluTable>1)]
						for(e in Mult){
							OtherAnnot=OtherAnnot[-which(OtherAnnot[,2]==e & OtherAnnot$as_type=="exclusion"),]
						}
					}
					if(any(table(OtherAnnot[,1])>1)){
						Mult=names(table(OtherAnnot[,1]))[which(table(OtherAnnot[,1])>1)]
						for(e in Mult){
							while(nrow(OtherAnnot[which(OtherAnnot[,1]==e),])>1){
								Del=which(OtherAnnot[,1]==e & OtherAnnot$as_type=="exclusion")
								OtherAnnot=OtherAnnot[-Del[1],]
							}
						}
					}
					
					ConvertNum=as.integer(substr(OtherAnnot[,1],6,nchar(OtherAnnot[,1])))-as.integer(substr(as.character(ASPSR),6,nchar(as.character(ASPSR))))
					names(ConvertNum)=OtherAnnot[,1]
					
					Left=ConvertNum[which(ConvertNum<0)]
					ClosestLeft=names(Left)[which.min(abs(Left))]
					
					Right=ConvertNum[which(ConvertNum>0)]
					ClosestRight=names(Right)[which.min(Right)]
					
					
					OtherAnnot=OtherAnnot[order(as.numeric(OtherAnnot$as_type)),]
					Junctions=OtherAnnot[OtherAnnot$PSR_ID%in%c(ClosestLeft,ClosestRight),2]
					if(length(Junctions)==2){
						if(length(unique(Annot[which(Annot[,1]%in%Junctions),2]))==1){
							Junctions=Junctions[1]
						}
					}
#					if(length(ClosestRight)==0){
#						if(is.null(JUCLeft)& !(is.null(JUCRight))){
#							Junctions=c(Junctions,JUCRight)
#						}
#						else if (!(is.null(JUCLeft)) & is.null(JUCRight)){
#							Junctions=c(Junctions,JUCLeft)
#						}
#					}
#					else if(length(ClosestLeft)==0){
#						if(is.null(JUCLeft) & !(is.null(JUCRight))){					
#							Junctions=c(Junctions,JUCRight)
#						}
#						
#						else if (!(is.null(JUCLeft)) & is.null(JUCRight)){
#							Junctions=c(Junctions,JUCLeft)
#						}
#					}
					#print(Junctions)
					
					if(length(Junctions)==0){
						Juc=Annot2[which(Annot2$as_type%in%c(3,5)),1][1]
						Junctions=c(Juc)
					}
					
					Del=c()
					for(j in Junctions){
						State=AnnotData[which(AnnotData$JUC_ID==j),ncol(AnnotData)][1]
#						if(State=="filtered out by INI"){
#							Del=c(Del,j)
#						}
					}
					if(length(Del)>0){
						if(length(Del)==length(Junctions)){
							Row=c(as.character(ASPSR),Junctions[1],Junctions[2],100000,100000,"INI Filtered Junctions")
							return(Row)
						}
						else{
							Row=c(as.character(ASPSR),Junctions[1],Junctions[2],100000,100000,"INI Filtered Junction")
							return(Row)
							
						}
					}
					
					
					#}
					
				}
				else{
					Annot=Annot[order(as.numeric(Annot$as_type)),]
					Junctions=Annot$JUC_ID[which(Annot$as_type%in%c("3","5"))]
					Del=c()
					for(j in Junctions){
						State=AnnotData[which(AnnotData$JUC_ID==j),ncol(AnnotData)][1]
#						if(State=="filtered out by INI"){
#							Del=c(Del,j)
#						}
					}
					if(length(Del)>0){
						if(length(Del)==length(Junctions)){
							Row=c(as.character(ASPSR),Junctions[1],Junctions[2],100000,100000,"INI Filtered Junctions")
							return(Row)
						}
						else{
							Junctions=Junctions[-which(Junctions%in%Del)]
							Row=c(as.character(ASPSR),Junctions[1],Junctions[2],100000,100000,"INI Filtered Junction")
							return(Row)
						}
					}
				}
				
				Juc=Annot$JUC_ID[which(Annot$as_type%in%c("3","5")&Annot$JUC_ID%in%Junctions)]
				
				PSR=as.vector(as.matrix(Data[which(Data$ExonID==as.character(ASPSR)),-c(1,2)]))
				PSR_Ranks=sort(PSR,index.return=TRUE)$ix
				JUC=as.vector(as.matrix(Data[which(Data$ExonID==Juc),-c(1,2)]))
				JUC_Ranks=sort(JUC,index.return=TRUE)$ix
				
				
				if(!Juc%in%unique(Data$ExonID)){
					Row=c(as.character(ASPSR),Juc,"",100000,100000,"Junction not found in Data")
					return(Row)
				}
				
				L=list(PSR_Ranks,JUC_Ranks)
				if(length(unique(c(length(PSR),length(JUC))))>1){
					MinLength=min(c(length(PSR),length(JUC)))
					Index=which(c(length(PSR),length(JUC))>MinLength)
					for(s in Index){
						L[[s]]=L[[s]][-which(L[[s]]%in%c((MinLength+1):length(L[[s]])))]
					}
				}
				
				Ranks=cbind(L[[1]],L[[2]])
				
				Y=as.vector(as.matrix(Ranks))
				exon=c(rep(1,nrow(Ranks)),rep(2,nrow(Ranks)))
				tissue=rep(c(rep(1,nrow(Ranks)/2),rep(2,nrow(Ranks)/2)),2)
				
				exon=as.factor(exon)
				tissue=as.factor(tissue)
				ft1<-lm(Y~exon*tissue-1)
				ft2<-lm(Y~exon+tissue)
				
				Chisq=lrtest(ft2,ft1)[2,4]
				
				InteractionSign=anova(ft1)[3,5]
				
				Row=c(as.character(ASPSR),Juc,"",InteractionSign,Chisq,"Only 1 type of junction")
				return(Row)
			}
		}
		else{
			
			if(all(!is.na(Annot$as_type))){
				if(all(Annot$as_type=="exclusion")){
					Row=c(as.character(ASPSR),"","",100000,100000,"Only exclusion junctions")
					return(Row)
				}
			}
			else{		
				Row=c(as.character(ASPSR),"","",100000,100000,"INI Filtered Junctions")
				return(Row)
			}
			
		}
		return(Row)
	}
	
	ScoresOutput=t(sapply(1:length(ASProbeSets), function(x) Scores(x,ASPSR=ASProbeSets[x],Annot=AnnotData[which(AnnotData$PSR_ID==ASProbeSets[x]),c(8,15)],AnnotData,Data)))
	return(ScoresOutput)
}


REIDS_RankedProbesets=ScoringBasedOnRanking(ASProbeSets=unique(ASPSR_PSR$ExonID), AnnotData=HTAData_RASA_PSRandJUN_LiverVSMuscle_AllPSR,Data=HTAData_RASA)
save(REIDS_RankedProbesets,file="REIDS_RankedProbesets.RData")




