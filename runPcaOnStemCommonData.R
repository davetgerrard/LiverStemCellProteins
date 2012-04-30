##runPcaOnStemCommonData.R


## requires loading of data with createStemCommonData.R

if(!exists('output')) { output <- FALSE } # would be over-ridden if already specified


excludeIndex <- which(rowMeans(stemcommondata[,dataColumns]) > 5)		# which rows are outliers.
cat("These rows are probable outliers\n")
print(stemcommondata[excludeIndex, ])

#dataColumns <- c(grep('H7',names(stemcommondata)),grep('H9',names(stemcommondata)) )
stem.pca <- prcomp(stemcommondata[,dataColumns])


if(output)  {

stem.pca <- prcomp(stemcommondata[,dataColumns])
pdf(file="plots/stemPCA.prcomp.pdf", width=10,height=10)
plot(stem.pca, main="Stem Cell Data\nVariance explained by Principal Components")
biplot(stem.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()

stem.pca <- prcomp(stemcommondata[,dataColumns], scale.=TRUE)
pdf(file="plots/stemPCA.prcomp.cor.pdf", width=10,height=10)
plot(stem.pca, main="Stem Cell Data\nVariance explained by Principal Components")
biplot(stem.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()

stem.pca <- princomp(stemcommondata[,dataColumns], cor=TRUE)		# this one gave a fantastic plot of PC2 vs PC3  (but is it biology or experimenatl effect?)
pdf(file="plots/stemPCA.princomp.scale.pdf", width=10,height=10)
plot(stem.pca, main="Stem Cell Data\nVariance explained by Principal Components")
biplot(stem.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()
	

stem.pca <- princomp(stemcommondata[-excludeIndex ,dataColumns], cor=TRUE)		
pdf(file="plots/stemPCA.princomp.cor.excludeOutliers.pdf", width=10,height=10)
plot(stem.pca, main="Stem Cell Data\nVariance explained by Principal Components")
biplot(stem.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()

stem.pca <- prcomp(stemcommondata[-excludeIndex ,dataColumns], scale.=TRUE)		
pdf(file="plots/stemPCA.prcomp.scale.excludeOutliers.pdf", width=10,height=10)
plot(stem.pca, main="Stem Cell Data\nVariance explained by Principal Components")
biplot(stem.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(stem.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=stemcommondata$Uniprot.Id[-excludeIndex]) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()


}	# end of if(output) block