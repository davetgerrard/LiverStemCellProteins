## Cliff Rowe's filtering script to select useable data columns from iTRAQ data files.
## Tweaked by DTG

H7rawdata<-read.csv(file="data/H7StemiTRAQ.csv") 

H7<-H7rawdata[c(2,8,5,6,9,14,19,24,29,34,39)]

H7pep2 <- H7[which(H7[,2]>1),]		### info: select data which was identified with 2 or more peptides
H795 <- H7pep2[which(H7pep2[,1]>1),]	### unused>1 (90% confidence in ID)

H7pep1 <- H7[which(H7[,2]==1),]		### select data which was identified with 1 peptide
H799 <- H7pep1[which(H7pep1[,1]>2),]	### 1 peptide (99% confidence in ID)
H7data<-rbind(H795[,3:11],H799[,3:11])	### dataset to use


names(H7data)<-c("Accession", "Name","H7.1","H7_HLC.1","H7.2","H7_HLC.2","H7.3","H7_HLC.3","H7_HLC.4")

H7data$Uniprot.Id <- matrix(unlist(strsplit(as.character(H7data$Accession),"|",fixed=T)),ncol=3,byrow=T)[,3]	# added by DTG to retain human readable unprot IDs
H7data$Uniprot.Accession <- matrix(unlist(strsplit(as.character(H7data$Accession),"|",fixed=T)),ncol=3,byrow=T)[,2]

write.csv(H7data,file="data/H7data.csv", row.names=F)
############################################


H9rawdata<-read.csv(file="data/H9StemiTRAQ.csv") 

H9<-H9rawdata[c(2,8,5,6,9,14,19,24,29,34,39)]

H9pep2 <- H9[which(H9[,2]>1),]		### 2+peptides
H995 <- H9pep2[which(H9pep2[,1]>1),]	### unused>1 (90%)

H9pep1 <- H9[which(H9[,2]==1),]		### 1 peptide
H999 <- H9pep1[which(H9pep1[,1]>2),]	### 1 peptide (99%)
H9data<-rbind(H995[,3:11],H999[,3:11])  ### remove unsed and peptides columns


names(H9data)<-c("Accession", "Name","H9.1","H9_HLC.1","H9.2","H9_HLC.2","H9.3","H9_HLC.3","H9_HLC.4")

H9data$Uniprot.Id <- matrix(unlist(strsplit(as.character(H9data$Accession),"|",fixed=T)),ncol=3,byrow=T)[,3]	# added by DTG to retain human readable unprot IDs
H9data$Uniprot.Accession <- matrix(unlist(strsplit(as.character(H9data$Accession),"|",fixed=T)),ncol=3,byrow=T)[,2]

write.csv(H9data,file="data/H9data.csv", row.names=F)