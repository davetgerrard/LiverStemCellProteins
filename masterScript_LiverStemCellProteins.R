## This script controls the overall analysis pipeline.

output <- FALSE
#output <- TRUE	#set this to TRUE if you want output of plots.

setwd("C:/Users/Dave/HanleyGroup/LiverStemCellProteins")

library(topGO)
library(qvalue)
library(gplots)

#source("filter_iTraq.R" )		# Filtering provided by Cliff. Based on their MS criteria for accepting proteins as detected.
						#only needs to be run once.



source("createStemCommonData.R")	# also outputs pair plots
## need to create different version of data and/or indexes
#	all samples, differentiated only
# 	Would be good to perform GO detect test using the full list of proteins detected in the differentiated samples.
source("prepareDataSubsets.R")


#	source("filterOutliers.R")	# would this be separate to other indexing?


## Prepare for GO enrichment tests. 
source("loadProt2Go.R")
source("C:/Users/dave/LiverProteins/scripts/topGoLiverProteinsUtilities.R")
source("detectGoStemCommonData.R")		# current version stop before incorporating func clustering.

source("detectGoStemHlcSamplesCommon.R")

##
source("runPcaOnStemCommonData.R")	# also ouptuts biplots

source("scoreGoStemCommonData.R")	# outputs a table of GO score results for each of first 5 PCs.
#source("")
#source("")