#################################################################################################################
####### This function identifies differential splicing and non-reliable probes in Affymetrix GeneChips	   ######
#################################################################################################################

#### This function computes array and exon score. 

gds<- function(SubgeneData, nsim=100) {
	
	nprobes <- nrow(SubgeneData)
	nsamples<- ncol(SubgeneData)
	nexon <- G <- length(unique(rownames(SubgeneData)))
	totalObservations <- nprobes*nsamples
	geneExp<- as.vector(as.matrix(SubgeneData))
	probes <- as.factor(as.vector(matrix(c(1:nprobes),nrow=nprobes,ncol=nsamples)))
	grp <- as.vector(matrix(rownames(SubgeneData),nrow=nprobes,ncol=nsamples))
	id <-   as.factor(as.vector(t(matrix(c(1:nsamples),nrow=nsamples,ncol=nprobes))))
	
	if(length(levels(probes))==1){
		geneRes=rep(0,length(probes))
	}
	else{
		ft <- lm(geneExp~probes+id-1,x=TRUE)  #Difference from inigds: id(samples) is involved in the calculation
	
	beta<- as.vector(summary(ft)$coefficient[,1])
	fixedDesignMatrix <- as.matrix(as.data.frame(ft$x)) ## create designsamples matrix
	geneRes <- as.vector(summary(ft)$residuals)
	}
	geneResMat <- matrix(geneRes,nprobes,nsamples)
	zid <- matrix(c(1:nsamples),nrow=nsamples,ncol=G)
	uz <- unique(grp)
	zcls2<- sapply(uz,function(x)ifelse(grp==x,1,0))
	idv <- as.numeric(id)
	Z<-matrix(0,totalObservations ,G*nsamples)
	for (ii in 1:nsamples) Z[idv==ii,(G*(ii-1)+1):(G*ii)]<-zcls2[idv==ii,]
	
	###########
	# Priors  #
	###########
	d0<-g0<-.001			
	tau <- 1
	
	#################
	# Store Results #
	#################
	nused <- floor(nsim/2)
	nburnin <-  nsim-nused
	outputbik <- matrix(0,nused,nsamples*G)
	outputerror <- NULL
	outputsigma <- matrix(0,nused,G)	# Fixed Effects
	###############################
	# Fixed Posterior Hyperparameter #
	#    for tau and taub		#
	###############################
	d<-d0+totalObservations/2
	
	###################
	# GIBBS SAMPLER   #
	###################
	
	nu0 <- G
	c0<-diag(G)			# Scale matrix for Wishart Prior onsamples Sigma.b
	nu<-nu0+nsamples
	Taub<-diag(G)			# Ransamplesdom Effects Prec Matrix
	z <- grp
	Taub <- diag(1,G,G)
	
	for (i in 1:nsim){
		
		set.seed(123*i+i)
		
		# Update b
		cZZ <- diag(colSums(Z))
		#vb <- diag(nsamples)%x%Taub+tau*cZZ #This is Var =D^-1+ sigma^(-2)*ni)
		vb <- solve(diag(nsamples)%x%Taub+tau*cZZ) # This is Var^(-1) =(D^-1+ sigma^(-2)*ni)^-1
		#vbInv=solve(vb)
		mb2<-(tau*crossprod(Z,geneRes))
		mb <- apply(vb,1,function(x) sum(x*mb2))
		vb1 <- diag(vb)
		b<-sapply(c(1:(nsamples*G)),function(x)rnorm(1,mean=mb[x],sd=sqrt(vb1[x])))
		bmat<-matrix(b,ncol=G,byrow=T)  		# Put insamples nsamples x G matrix form for updatinsamplesg taub
		
		
		# Update tau
		zb <- rowSums(sapply(c(1:G),function(x) t(matrix(bmat[,x],nsamples,nprobes))*matrix(zcls2[,x],nprobes,nsamples)))	
		geneReszb <-geneRes-zb	
		g<-g0+crossprod(geneReszb,geneReszb)/2
		tau<-rgamma(1,d,g)
		
		
		# Update Taub 
		
		Sigma.b<-riwish(nu,c0+crossprod(bmat,bmat))
		Taub<-solve(Sigma.b)
		
		if(i > nburnin ){ 
			outputbik[(i-nburnin),] <- b
			outputerror[(i-nburnin)] <- 1/tau
			outputsigma[(i-nburnin),] <- as.vector(diag(Sigma.b))
			
		}
	}
	
	tp1 <- apply(outputsigma,2,function(x)quantile(x, probs = c(0.025, 0.5, 0.95)))
	tp2 <- quantile(outputerror, probs = c(0.025, 0.5, 0.95))
	icc <- t(apply(tp1,2,function(x) x/(x+tp2)))
	icc <- data.frame(exon=uz,type="icc",icc )
	nik <- as.vector(t(sapply(uz,function(x) rep(x,nsamples))))
	ikEffect2  <- t(apply(outputbik,2,function(x)quantile(x, probs = c(0.025, 0.5, 0.95))))
	ikEffect<- data.frame(exon=nik,type="bik",ikEffect2)
	exonscore <- icc[,c(1,4)]
	arrayscore <- ikEffect[,4] 
	
	dim(arrayscore)<- c(G,nsamples)
	rownames(arrayscore)<-ikEffect$exon[1:G] 
	colnames(arrayscore) <- colnames(SubgeneData)
	Output <- list(exonScores =exonscore ,arrayScores=arrayscore)
	return(Output)
}