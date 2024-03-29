

if(!exists('output')) { output <- FALSE }  # would be over-ridden if already specified

###################INFO: Re-run GO analyses grouping GO terms by functional clusters.

###INFO:  Go through func clusters and pull in topGO results into a single table across BP,MF,CC

# standard topGO settings
elimCutOff <- 0.01
#pcs <- 1:4
goGraphs <- c("BP","MF","CC")
nodeSizeValue <- 10	
topNodesValue <- 20
topTerms <- 50

###INFO: set up a topGOdata object for each ontology based on the list of interesting proteins -> goDataCollection.all
## Used in performing GO tests and to retrieve appropriate protein lists from GO terms. 


allProtList <- names(prot2go)
##INFO: in this GO analysis use all proteins detected in at least one sample. N.B. there are many with zero data.
detectedProts.Hlc.all <- stemcommon$Uniprot.Accession[row_index.HLC.all]
geneList.all.binary <- factor(as.integer(allProtList %in% detectedProts.Hlc.all ))
names(geneList.all.binary) <- allProtList

goDataCollection.all.binary <- list()
for(thisGOgraph in goGraphs )  {
	goDataCollection.all.binary[[thisGOgraph]] <- new("topGOdata",
		description =  "Stem cell Hlc Samples all proteins data set",
		ontology = thisGOgraph,
		allGenes = geneList.all.binary,
		nodeSize = nodeSizeValue ,
		annot = annFUN.GO2genes,
		GO2genes=go2prot 
		)
}




#INFO:Base GO analysis -> "liverProteinsUbiqBasicGo.all.pdf"

#pdf(file="liverProteinsUbiqBasicGo.all.pdf", width=10,height=9)
### perform the GO analysis, iterating through ontologies and binding all results into one table. 
summaryDetectResults.all <- data.frame()
for(thisGOgraph in goGraphs )  {
	summaryDetectResults.all <- rbind(summaryDetectResults.all ,runDetectGoTests(goDataCollection.all.binary[[thisGOgraph]]))
}
#dev.off()

######## END OF BASE GO ANALYSIS


#INFO: output base GO results ->
write.table(summaryDetectResults.all[order(summaryDetectResults.all$elimFisher),],file="output/detectGOSummary.HLC.all.tab",quote=F,row.names=F,sep="\t")




stopifnot(FALSE)	#  CURRENTLY STOPS HERE.	 The remainder used func clustering results.



### CALL to functional clustering. 
#source("C:/Users/dave/LiverProteins/scripts/funcClusterGoTable.R")


#INFO:  Join GO terms into a group if the functional cluster contains a high proportion of the  proteins linked to that GO term
goContainProp <- 0.8  # same result with 0.7
detectTableSigThreshold <- 1.0e-05

#INFO:  ?need to order the results (and start at the top?)
# perhaps filter by significance cut-off
orderBy <- "elimFisher"



#INFO: NB clusters are built only with ubiquitous proteins. Only going to use these proteins in clustering. 
validProts <- proteinByLiverSample.ubiquitous$spAccession

# find all relevant clusters for each GO term. Not needed
#summaryDetectResults.all$allClusters  <- apply(summaryDetectResults.all,1, FUN= function(x) listAllOverlappingClusters(
#						goProteins.valid=intersect(listProtsInGoFromList(goTerm=x["goTerm"],ontology=x["ontology"],goDataCollection.all.binary),validProts),#
#						seedList	,goContainProp ))

#INFO: find the best cluster for a GO term.
summaryDetectResults.all$bestCluster  <- apply(summaryDetectResults.all,1, FUN= function(x) listBestOverlappingCluster(
						goProteins.valid=intersect(listProtsInGoFromList(goTerm=x["goTerm"],ontology=x["ontology"],goDataCollection.all.binary),validProts),
						seedList	,goContainProp ))

#summaryDetectResults.all <- subset(summaryDetectResults.all, select=-c(allClusters, bestCluster))	# testing, needed to undo above.


write.table(summaryDetectResults.all[order(summaryDetectResults.all$elimFisher),],file="detectGOWithClusterNumberElimFisher.all.tab",sep="\t",quote=F,row.names=F)

#INFO: use clusterTable() to filter and cluster the table
summaryDetectResults.all.sigClustered <- clusterTable(summaryDetectResults.all,orderBy=orderBy,detectTableSigThreshold=detectTableSigThreshold)     # function(table,orderBy,detectTableSigThreshold = 1.0e-05)  {

write.table(summaryDetectResults.all.sigClustered,file="sigDetectGoTableByCluster.all.tab",quote=F,row.names=F,sep="\t")




