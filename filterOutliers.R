###filterOutliers.R


## after running PCA on each set separately, eye-balled following strong outliers:-
#H7	TIM50_HUMAN, CO2A1_HUMAN
#H9	RIC3_HUMAN, MILK2_HUMAN, CRYAA_HUMAN


## need to remove these from dataset before rest of analyses.

#inspect the relevant results

stemcommon[match(c("TIM50_HUMAN", "CO2A1_HUMAN","RIC3_HUMAN", "MILK2_HUMAN", "CRYAA_HUMAN"),stemcommon$Uniprot.Id),]


maxAboveTen <- apply(stemcommon[,dataColumns], 1, max, na.rm=T) > 10