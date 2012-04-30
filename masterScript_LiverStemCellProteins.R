## This script controls the overall analysis pipeline.

output <- FALSE
#output <- TRUE	#set this to TRUE if you want output of plots.

setwd("C:/Users/Dave/HanleyGroup/LiverStemCellProteins")



#source("filter_iTraq.R" )		# only needs to be run once.


source("createStemCommonData.R")

source("runPcaOnStemCommonData.R")

#source("")