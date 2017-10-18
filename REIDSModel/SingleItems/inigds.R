inigds<- function(SubgeneData, nsim=1000) {
	
	# Preparation 
	nprobes <- nrow(SubgeneData)
	nsamples<- ncol(SubgeneData)
	nexon <- G <- length(unique(rownames(SubgeneData)))
	totalObs <- nprobes*nsamples
	geneExp<- as.vector(as.matrix(SubgeneData))
	probes <- as.factor(as.vector(matrix(c(1:nprobes),nrow=nprobes,ncol=nsamples)))
	grp <- as.vector(matrix(rownames(SubgeneData),nrow=nprobes,ncol=nsamples))
	id <-   as.factor(as.vector(t(matrix(c(1:nsamples),nrow=nsamples,ncol=nprobes))))
	
	# Fit of model
	if(length(levels(probes))==1){
		geneRes=rep(0,length(probes))
	}
	else{
		ft <- lm(geneExp~probes-1,x=TRUE) # Why is the model fitted? To get resiudal for estimation of random intercepts?
	
		beta<- as.vector(summary(ft)$coefficient[,1]) #not used?
		fixedDesignMatrix <- as.matrix(as.data.frame(ft$x)) ## create designsamples matrix: not used?
		geneRes <- as.vector(summary(ft)$residuals)  #get residuals ==> for estimation of the ramdom intercepts??
	}
	geneResMat <- matrix(geneRes,nprobes,nsamples) #filled in by column: ok, one column per sample, calculation happens per sample
	
	zid <- matrix(c(1:nsamples),nrow=nsamples,ncol=G)
	uz <- unique(grp)
	zcls2<- sapply(uz,function(x)ifelse(grp==x,1,0))
	
	idv <- as.numeric(id)
	Z<-matrix(0,totalObs ,G*nsamples)
	for (ii in 1:nsamples) Z[idv==ii,(G*(ii-1)+1):(G*ii)]<-zcls2[idv==ii,]
	
	###########
	# Priors  #
	###########
	d0<-g0<-.001		# alpha and beta	
	tau <- 1 #sigma^-2 of measurement error
	
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
	d<-d0+totalObs/2  #alpah+0.5*N
	
	###################
	# GIBBS SAMPLER   #
	###################
	
	nu0 <- G
	c0<-diag(G)			# Scale matrix for Wishart Prior onsamples Sigma.b
	nu<-nu0+nsamples   #gamma+n 
	Taub<-diag(G)			# Ransamplesdom Effects Prec Matrix
	z <- grp
	Taub <- diag(1,G,G) #not needed
	
	for (i in 1:nsim){
		
		set.seed(123*i+i)
		
		# Update b
		cZZ <- diag(colSums(Z))
		#vb <- diag(nsamples)%x%Taub+tau*cZZ #This is Var =D^-1+ sigma^(-2)*ni)
		vb <- solve(diag(nsamples)%x%Taub+tau*cZZ) # This is Var^(-1) =(D^-1+ sigma^(-2)*ni)^-1
		#vbInv=solve(vb)
		mb2<-(tau*crossprod(Z,geneRes ))  #Theta
		mb <- apply(vb,1,function(x) sum(x*mb2) ) #Var^-1*Theta
		vb1 <- diag(vb)  #only need variances here
		b<-sapply(c(1:(nsamples*G)),function(x)rnorm(1,mean=mb[x],sd=sqrt(vb1[x])))  #aren't these the variances^-1 instead of the variances?
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
	
	tp1 <- apply(outputsigma,2,function(x) mean(x))
	tp2 <- mean(outputerror)
	icc <- tp1/(tp1+tp2)
	
	Output <- icc
	return(Output)
}