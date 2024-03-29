

########### PRE-REQUISITES



# A map of protein to GO annotations. 

# for the funcitonal clustering, the proteins need to have been clustered into functional groups based on shared annotations.

# load utilities for performing functional clustering of GO result tables.

# load utilities for running multiple GO analyses.



########### PROCESS


if(!exists('output')) { output <- FALSE }  # would be over-ridden if already specified

###################INFO: Re-run GO analyses grouping GO terms by functional clusters.

###INFO:  Go through func clusters and pull in topGO results into a single table across BP,MF,CC

# standard topGO settings
elimCutOff <- 0.01
pcs <- 1:4
goGraphs <- c("BP","MF","CC")
nodeSizeValue <- 10	
topNodesValue <- 20
topTerms <- 50

#INFO: set up a topGOdata object for each ontology based on the list of interesting proteins -> goDataCollection.ubiq
#INFO: Used in performing GO tests and to retrieve appropriate protein lists from GO terms. 


stem.princomp.cor.scores <- cbind(subset(stemcommondata, select=c("Name","Uniprot.Id", "Uniprot.Accession")) ,stem.princomp.cor$scores)


geneList.ubiq <- stem.princomp.cor.scores$Comp.1			#take an aribitrary set of scores.
names(geneList.ubiq) <- stem.princomp.cor.scores$Uniprot.Accession	# the names are important.

goDataCollection.ubiq.scored <- list()
for(thisGOgraph in goGraphs )  {
	goDataCollection.ubiq.scored[[thisGOgraph]] <- new("topGOdata",
		description =  "Stem Cell proteins data set",
		ontology = thisGOgraph,
		allGenes = geneList.ubiq,
  		geneSelectionFun = topDiffGenes,
		nodeSize = nodeSizeValue ,
		annot = annFUN.GO2genes,
		GO2genes=go2prot 
		)
}

####INFO: Perform tests and build set of result tables. One per PC.
#INFO: following structure iterates through PCs and within each through GO ontologies (BP, MF, CC)
#INFO: runs multiple tests for GO enrichment of PC scores and stores p-values per GO term
#INFO: aim is to create one table for each PC with the results from the ontologies combined. 

summaryScoreResults <- list()	#INFO: the PCs are treated separately. Their individual results are stored in elements of a list. 

for( i in 1:5)  {
	compHead <- paste("Comp.",i,sep="")		# this PC
	geneList <- stem.princomp.cor.scores[,compHead]	# get PC scores for proteins
	names(geneList) <- stem.princomp.cor.scores$Uniprot.Accession		# make sure scores have names of proteins

	summaryScoreResults[[i]] <- data.frame()	#initiate the data frame for all results from this PC

	for(thisGOgraph in goGraphs )  {		# loop through 3 GO ontologies
		#INFO: The same goData object is re-used but must be updated with the relevant scores before use.
		goDataCollection.ubiq.scored[[thisGOgraph]] <- updateGenes(goDataCollection.ubiq.scored[[thisGOgraph]],geneList,topDiffGenes)		
		#INFO: runScoreGoTests() performs the tests and returns an ouptut table
		summaryScoreResults[[i]] <- rbind(summaryScoreResults[[i]],runScoreGoTests(goDataCollection.ubiq.scored[[thisGOgraph]],geneList,topDiffGenes))

	}
}

#INFO: output the plain GO on PC results
for(i in 1:length(summaryScoreResults)) {
	outFileName <- paste("GoSummaryByPc",i,"tab", sep=".")
	write.table(summaryScoreResults[[i]],file=paste("output",outFileName,sep="/"),sep="\t",quote=F,row.names=F) 
	rm(outFileName)
}


stopifnot(FALSE)


#INFO:  Join GO terms into a group if the functional cluster contains a high proportion of the  proteins linked to that GO term
goContainProp <- 0.8  # same result with 0.7
detectTableSigThreshold <- 1.0e-03

#INFO:  ?need to order the results (and start at the top?)
# perhaps filter by significance cut-off
orderBy <- "elimWilcox"
#INFO: Assign best functional cluster to GO terms, if applicable.
###PC table...

#summmaryPcResultList[[i]]

#INFO: output the GO results with clusters marked on PC results
for(i in 1:length(summaryScoreResults)) {
	
	summaryScoreResults[[i]]$bestCluster  <- apply(summaryScoreResults[[i]],1, FUN= function(x) listBestOverlappingCluster(
						goProteins.valid=intersect(listProtsInGoFromList(goTerm=x["goTerm"],ontology=x["ontology"],goDataCollection.ubiq.scored),validProts),
						seedList	,goContainProp ))
	outFileName <- paste("scoreGoWithClusterByPc",i,"tab", sep=".")
	write.table(summaryScoreResults[[i]][order(summaryScoreResults[[i]][,orderBy]),],file=outFileName,sep="\t",quote=F,row.names=F) 	
	rm(outFileName)
}


# use same detectTableSigThreshold for each PC table?


#INFO: Use clusterTable() to ouput a filtered and clustered summary of the top terms in each PC table.
summaryScoreResults.sigClustered <- list()
for(i in 1:length(summaryScoreResults)) {
	summaryScoreResults.sigClustered[[i]] <- clusterTable(summaryScoreResults[[i]],orderBy=orderBy,detectTableSigThreshold=detectTableSigThreshold)     # function(table,orderBy,detectTableSigThreshold = 1.0e-05)  {
	outFileName <- paste("GoSummaryTopClusteredByPc",i,"tab", sep=".")
	write.table(summaryScoreResults.sigClustered[[i]],file=outFileName,sep="\t",quote=F,row.names=F) 	
	rm(outFileName)
}




