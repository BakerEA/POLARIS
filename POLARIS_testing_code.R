
library(data.table)

for (sim in 1:10){


	eff<- fread(paste("./POLARIS_test_data/ScenarioA1/Discovery/simulation", sim, ".assoc.logistic", sep=""), header=T, data.table=FALSE)
	
	dat<- fread(paste("./POLARIS_test_data/ScenarioA1/Test/simulation", sim, ".raw", sep=""), header=T, data.table=FALSE)

    #Remove last two digits from header to match summ file
    header<- substr(colnames(dat)[7:ncol(dat)], 1, nchar(colnames(dat)[7:ncol(dat)])-2)
    colnames(dat)[7:ncol(dat)]<- header
    
    #Subset SNPs in the summary dataset
    a<- which(colnames(dat) %in% eff$SNP)
    dat<- dat[, c(1:6,a)]
	
	
	for (i in 7:ncol(dat)){
		
		s<- which(eff$SNP==colnames(dat)[i])
		
		#Change the 9 to whichever column is MAF
		dat[is.na(dat[,i]),i]<- 2*eff[s,9]
		
	}
	 
	co<-cor(as.matrix(dat[,7:ncol(dat)]))
 
	ee<-eigen(co)
	la0<-1/sqrt(nrow(dat))
 
	PM<-ee$vectors %*% diag(sqrt((1+la0)/(ee$values+la0))) %*% t(ee$vectors)
	POLARIS<- eff$BETA %*% PM %*% t(data.matrix(dat[,7:ncol(dat)]))
					
	fit<- glm((dat$PHENOTYPE-1)~t(POLARIS), family=binomial)
	summary(fit)

}