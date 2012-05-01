### A place to test things


nrow(H7data)

nrow(H9data)

nrow(stemcommon) # had all proteins detected.


# restrict to culture samples.

col_index.HLC <- grep('HLC', names(stemcommon))

stemcommon.HLC <- cbind(subset(stemcommon, select=c("Accession","Name","Uniprot.Id","Uniprot.Accession")), stemcommon[,col_index.HLC])

HLC.allNA <- which(rowSums(is.na(stemcommon[,col_index.HLC])) != ncol(stemcommon[,col_index.HLC]))

row_index.HLC.all <- rowSums(is.na(stemcommon[,col_index.HLC])) != ncol(stemcommon[,col_index.HLC])	# which proteins are present in some of the HLC samples
row_index.HLC.common <- which(complete.cases(stemcommon[,col_index.HLC]))


## indexes
col_index.H9 <- grep('H9', names(stemcommon))
col_index.H7 <- grep('H7', names(stemcommon))
col_index.HLC <- grep('HLC', names(stemcommon))
col_index.ESC <- setdiff(union(col_index.H9,col_index.H7),col_index.HLC)

stemcommon[1,col_index.ESC]



H7.complete <- H7data[complete.cases(H7data),]
nrow(H7data)
nrow(H7.complete)
H7.pca <- princomp(H7.complete[,grep('H7',names(H7.complete))], cor=T)
H7.pca.cor.scores <- cbind(subset(H7.complete, select=c("Name","Uniprot.Id", "Uniprot.Accession")) ,H7.pca$scores)
write.table(H7.pca.cor.scores, file=paste("output","H7.pca.cor.scores.tab",sep="/"),row.names=F, quote=F, sep="\t")
pdf(file="plots/H7.PCA.princomp.cor.pdf", width=10,height=10)
plot(H7.pca, main="H7 Stem Cell Data\nVariance explained by Principal Components")
biplot(H7.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.pca,choices=c(1,6), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()

H7.complete.noOutliers <- H7.complete[apply(H7.complete[,grep('H7',names(H7.complete))], 1, max, na.rm=T) < 10, ]
nrow(H7.complete.noOutliers)
H7.noOutliers.pca <- princomp(H7.complete.noOutliers[,grep('H7',names(H7.complete.noOutliers))], cor=T)
H7.noOutliers.pca.cor.scores <- cbind(subset(H7.complete.noOutliers, select=c("Name","Uniprot.Id", "Uniprot.Accession")) ,H7.noOutliers.pca$scores)
write.table(H7.noOutliers.pca.cor.scores, file=paste("output","H7.noOutliers.pca.cor.scores.tab",sep="/"),row.names=F, quote=F, sep="\t")
pdf(file="plots/H7.noOutliers.PCA.princomp.cor.pdf", width=10,height=10)
plot(H7.noOutliers.pca, main="H7 Stem Cell Data\nVariance explained by Principal Components")
biplot(H7.noOutliers.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.noOutliers.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.noOutliers.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.noOutliers.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.noOutliers.pca,choices=c(1,6), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H7.noOutliers.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=H7.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()


H9.complete <- H9data[complete.cases(H9data),]
nrow(H9data)
nrow(H9.complete)
H9.pca <- princomp(H9.complete[,grep('H9',names(H9.complete))])
H9.pca.cor.scores <- cbind(subset(H9.complete, select=c("Name","Uniprot.Id", "Uniprot.Accession")) ,H9.pca$scores)
write.table(H9.pca.cor.scores, file=paste("output","H9.pca.cor.scores.tab",sep="/"),row.names=F, quote=F, sep="\t")
pdf(file="plots/H9.PCA.princomp.cor.pdf", width=10,height=10)
plot(H9.pca, main="H9 Stem Cell Data\nVariance explained by Principal Components")
biplot(H9.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.pca,choices=c(1,6), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()

H9.complete.noOutliers <- H9.complete[apply(H9.complete[,grep('H9',names(H9.complete))], 1, max, na.rm=T) < 10, ]
nrow(H9.complete.noOutliers)
H9.noOutliers.pca <- princomp(H9.complete.noOutliers[,grep('H9',names(H9.complete.noOutliers))], cor=T)
H9.noOutliers.pca.cor.scores <- cbind(subset(H9.complete.noOutliers, select=c("Name","Uniprot.Id", "Uniprot.Accession")) ,H9.noOutliers.pca$scores)
write.table(H9.noOutliers.pca.cor.scores, file=paste("output","H9.noOutliers.pca.cor.scores.tab",sep="/"),row.names=F, quote=F, sep="\t")
pdf(file="plots/H9.noOutliers.PCA.princomp.cor.pdf", width=10,height=10)
plot(H9.noOutliers.pca, main="H9 Stem Cell Data\nVariance explained by Principal Components")
biplot(H9.noOutliers.pca,choices=c(1,2), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.noOutliers.pca,choices=c(1,3), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.noOutliers.pca,choices=c(1,4), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.noOutliers.pca,choices=c(1,5), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.noOutliers.pca,choices=c(1,6), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
biplot(H9.noOutliers.pca,choices=c(2,3), col=c("grey","black"),cex=c(0.5,1),xlabs=H9.complete.noOutliers$Uniprot.Id) ;abline(v=0,lty=2) ; abline(h=0,lty=2) 
dev.off()

