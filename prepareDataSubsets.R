#prepareDataSubsets.R


## indexes
col_index.H9 <- grep('H9', names(stemcommon))
col_index.H7 <- grep('H7', names(stemcommon))
col_index.HLC <- grep('HLC', names(stemcommon))
col_index.ESC <- setdiff(union(col_index.H9,col_index.H7),col_index.HLC)

row_index.HLC.all <- rowSums(is.na(stemcommon[,col_index.HLC])) != ncol(stemcommon[,col_index.HLC])	# which proteins are present in some of the HLC samples
# check if these contain rows that are ALL NA.
row_index.HLC.common <- which(complete.cases(stemcommon[,col_index.HLC]))
