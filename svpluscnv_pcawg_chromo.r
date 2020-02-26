### Evaluation of performance of shattered.regions against PCAWG chromothripsis survey

# data:
# PCAWG all tumor types data contains CNVs and SVs derived from WGS
# download from GDC legacy portal (https://portal.gdc.cancer.gov/legacy-archive/
load("data/pcawg.rda",verbose=TRUE)

# data:
# Chromothripsis survey from PCAWG table contains info about complex chromosomal rearrangements in PCAWG samples
# download from manuscript (https://doi.org/10.1038/s41588-019-0576-7) based on ShatterSeek algorithm (https://github.com/parklab/ShatterSeek)
# suppl table 1 (https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0576-7/MediaObjects/41588_2019_576_MOESM3_ESM.xlsx)

load("data/pcawg_chromo.rda",verbose=TRUE)
SSeek_results <- pcawg_chromo[which(pcawg_chromo$chromo_label %in% c("High confidence","Low confidence") ),]


# run shattered.regions for all PCAWG samples
shreg_all_pcawg <- shattered.regions(cnv_pcawg, svc_pcawg, fc.pct = 0.1, interleaved.cut = 0.33, verbose=TRUE)

# collect shattered regions from the output
shreg_all_pcawg <- shattered.regions(cnv_pcawg, svc_pcawg, fc.pct = 0.1, interleaved.cut = 0.33, verbose=TRUE)
# in order to link the sampleids we are going to use donor_unique_id field
spcid <- unname(unlist(sapply(names(shreg_all_pcawg@regions.summary), 
                              function(i) rep(i, nrow(shreg_all_pcawg@regions.summary[[i]])) )))
duid <- ali2dui[spcid]  
SReg_results <- data.table(spcid,duid,do.call(rbind,shreg_all_pcawg@regions.summary))

# to make data comparable we restrict the analysis to samples included by the PCAWG study only (see https://doi.org/10.1038/s41588-019-0576-7 for filtering details) 
SReg_results_in <- SReg_results[which(SReg_results$duid %in% unique(pcawg_chromo$donor_unique_id)),]


# paste donor_id and chromosome from shattered regions from shattered.regions and PCAWG study
SSeek_chr <-  unique(gsub("_","_chr",unique(unite(SSeek_results, newcol, c(donor_unique_id,Chr), remove=FALSE,sep="_")$newcol)))
SReg_chr <-  unique(unite(SReg_results_in, newcol, c(duid,chrom), remove=FALSE,sep="_")$newcol)

#The total number of tests is defined by the number of samples in the gold standard times chromosomes (23):
tot <- length(unique(pcawg_chromo$donor_unique_id))*23
TP <- length(intersect(SSeek_chr,SReg_chr))
FP <- length(setdiff(SReg_chr,SSeek_chr))
FN <- length(setdiff(SSeek_chr,SReg_chr))
TN <- tot - TP - FP -FN

# measure different scores
Precision <- TP/(TP+FP)
Recall <- TP/(TP+FN)
F1 <- 2*Precision*Recall/(Precision+Recall)
MCC <- (TP*TN -FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
Sens <- TP/(TP+FN)
Spec <- TN/(TN+FN)
Accuracy = (TP+TN)/(TP+TN+FP+FN)


message(paste("Accuracy=",Accuracy,"\nPrecission=",Precision,"\nRecall=",Recall))
