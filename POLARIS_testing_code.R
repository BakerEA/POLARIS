
strt_time<-Sys.time()

countp05<-0
countp01<-0
countp001<-0

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
	
	#Include PRS
	
	prs<- snps %*% summ$BETA
	prs_norm<- ((prs-mean(prs))/sd(prs))
	
	
    data1<-data.frame((repl$PHENOTYPE-1), prs_norm)
    colnames(data1)[1]<- "Status"
	
	
    null1<- glm(Status~1, data=data1, family="binomial")
	fit1<- glm(Status~ prs_norm, data=data1, family="binomial")
	
	p1<-1-pchisq(null1$deviance-fit1$deviance,1)
	
	if (p1<=0.05){count1p05<-count1p05+1}
	if (p1<=0.01){count1p01<-count1p01+1}
	if (p1<=0.001){count1p001<-count1p001+1}
	
}


print("POLARIS")
print( "P<0.05")
print(countp05)
print( "P<0.01")
print(countp01)
print( "P<0.001")
print(countp001)

print(Sys.time()-strt_time)





