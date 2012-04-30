#createStemCommonData.R
if(!exists('output')) { output <- FALSE } # would be over-ridden if already specified

H7data <- read.csv(file="data/H7data.csv")
H9data <- read.csv(file="data/H9data.csv")


stemcommon<-merge(H7data,H9data,incomparables=NA,all=TRUE)	# use the default 'by' to prevent duplicate columns
stemcommondata<-stemcommon[complete.cases(stemcommon), ]

nrow(stemcommon)
nrow(stemcommondata)		# 1922 as per my method above.

dataColumns <- c(grep('H7',names(stemcommondata)),grep('H9',names(stemcommondata)) )
stemcommondata[1,dataColumns]


if(output)  {	
pdf("plots/stemCommonData.pairPlots.pdf")
pairs(stemcommondata[ ,dataColumns])
dev.off()

}