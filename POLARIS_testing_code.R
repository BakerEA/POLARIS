
strt_time<-Sys.time()

p_thr<- 1

for (sim in 1:10){

    repl<- read.table(paste("./POLARIS_test_data/ScenarioA1/N10000_OR11_20/Test/simulation", sim, ".raw", sep=""), header=T)

    summ<- read.table(paste("./POLARIS_test_data/ScenarioA1/N10000_OR11_20/Discovery/simulation", sim, ".assoc.logistic", sep=""), header=T)
	
    #Filter out inclusion SNPs
	summ<- summ[summ$P<=p_thr,]

    #Remove last two digits from header to match summ file
    header<- substr(colnames(repl)[7:ncol(repl)], 1, nchar(colnames(repl)[7:ncol(repl)])-2)
    colnames(repl)[7:ncol(repl)]<- header
    
    #Subset SNPs in the summary dataset
    a<- which(colnames(repl) %in% summ$SNP)
    repl<- repl[, c(1:6,a)]

	#Separate SNPs only
    snps<- as.matrix(repl[, 7:ncol(repl)])
	
	#Spectral decomposition of the correlation matrix
	spect_decomp<- eigen(cor(snps))
	
	eigval<- spect_decomp$values
	eigvec<- spect_decomp$vectors
	
	A<-  summ$BETA %*% eigvec
	
	eigwgt<- eigvec
	
	for (i in 1:ncol(snps)){
        eigwgt[,i]<- eigwgt[,i]* sqrt((1+(1/sqrt(nrow(snps))))/(eigval[i]+(1/sqrt(nrow(snps)))))
	}
	
    #Find POLARS adjusted Betas
	Badj<- eigwgt %*% t(A)
	
    #Compute POLARIS
	polaris<- snps  %*%  Badj
	
    #Normalise POLARIS scores
    polaris_norm<-((polaris-mean(polaris))/sd(polaris))
    
    #Create data with phenotype and normalised POLARIS scores
    data<-data.frame((repl$PHENOTYPE-1), polaris_norm)
    colnames(data)[1]<- "Status"
    
    #Fit model for score to find power
	null<- glm(Status~1, data=data, family="binomial")
	fit<- glm(Status~ polaris_norm, data=data, family="binomial")
	
	p<-1-pchisq(null$deviance-fit$deviance,1)
	
}

print(Sys.time()-strt_time)





