## This script controls the overall analysis pipeline.

output <- FALSE
#output <- TRUE	#set this to TRUE if you want output of plots.

setwd("C:/Users/Dave/HanleyGroup/LiverStemCellProteins")

library(topGO)
library(qvalue)
library(gplots)

#source("filter_iTraq.R" )		# only needs to be run once.


source("createStemCommonData.R")	# also outputs pair plots
## need to create different version of data
#	all samples, differentiated only
# 	Would be good to perform GO detect test using the full list of proteins detected in the differentiated samples.



#	source("filterOutliers.R")

source("loadProt2Go.R")
source("C:/Users/dave/LiverProteins/scripts/topGoLiverProteinsUtilities.R")
#	source("detectGoStemCommonData.R")

source("runPcaOnStemCommonData.R")	# also ouptuts biplots

#	source("pcScoreGoStemCommonData.R")
#source("")
#source("")