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

