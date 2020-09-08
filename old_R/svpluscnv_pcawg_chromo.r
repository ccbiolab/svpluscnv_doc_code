library(svpluscnv)
library(tidyr)
library(data.table)
library(VennDiagram)

# local dir
setwd("<path-to-local-folder>/svpluscnv_doc_code/")
setwd("~/Box Sync/git/svpluscnv_doc_code/")

# data:
# PCAWG all tumor types data contains CNVs and SVs derived from WGS
# download from GDC legacy portal (https://portal.gdc.cancer.gov/legacy-archive/
load("data/pcawg.rda",verbose=TRUE)


## CHROMOTHRIPSIS analysis of chromothripsis based on ShatterProof, ShatterSeek and shattered.regions (svpluscnv);
# comparison across methods

# load data from 
# Generated using code available at ~/shatterproof/prepare_shatterproof.r and ~/shatterproof/run_shatterproof.pl
load("data/shatterproof_pcawg_chromo_results.rda",verbose=TRUE)

## chromothripsis data from PCAWG
# Chromothripsis survey from PCAWG table contains info about complex chromosomal rearrangements in PCAWG samples
# download from manuscript (https://doi.org/10.1038/s41588-019-0576-7) based on ShatterSeek algorithm (https://github.com/parklab/ShatterSeek)
# suppl table 1 (https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0576-7/MediaObjects/41588_2019_576_MOESM3_ESM.xlsx)

load("data/pcawg_chromo.rda",verbose=TRUE)
SSeek_results <- pcawg_chromo[which(pcawg_chromo$chromo_label %in% c("High confidence","Low confidence") ),]
SSeek_results_breast <- SSeek_results[which(SSeek_results$histology_abbreviation == "Breast-AdenoCA" ),]

intersect(unique(do.call(rbind,strsplit(SSeek_results$donor_unique_id,"::"))[,2]),
names(ali2dui) )

sseek <- SSeek_results
SSeek.bed <- as.data.frame(SSeek_results[,c("donor_unique_id","Chr","Start","End")])
SSeek.bed[,"Chr"] <- paste("chr",SSeek.bed[,"Chr"] ,sep="")
colnames(SSeek.bed) <- c("sampleid","chrom","start","end")

# run shattered.regions for all PCAWG samples
#shreg_all_pcawg <- shattered.regions(cnv_pcawg, svc_pcawg, fc.pct = 0.1, interleaved.cut = 0.33, verbose=TRUE)
#save(shreg_all_pcawg, file="data/shreg_all_pcawg.rda")
load("data/shreg_all_pcawg.rda")
sregs <- shreg_all_pcawg

# in order to link the sampleids we are going to use donor_unique_id field
spcid <- unname(unlist(sapply(names(sregs@regions.summary), 
                              function(i) rep(i, nrow(sregs@regions.summary[[i]])))))
duid <- ali2dui[spcid]  
SReg_results <- data.table(spcid,duid,do.call(rbind,sregs@regions.summary))

# to make data comparable we restrict the analysis to samples included by the PCAWG study only (see https://doi.org/10.1038/s41588-019-0576-7 for filtering details) 
SReg_results_in <- SReg_results[which(SReg_results$duid %in% unique(pcawg_chromo$donor_unique_id)),]
dataf <- data.frame(unlist(unname(sapply(names(shreg_all_pcawg@regions.summary), function(i) rep(i,nrow(shreg_all_pcawg@regions.summary[[i]]))))),
                    do.call(rbind,shreg_all_pcawg@regions.summary))
sregs.bed <- dataf[,1:4]
colnames(sregs.bed) <- c("sampleid","chrom","start","end")

# paste donor_id and chromosome from shattered regions from shattered.regions and PCAWG study
SSeek_chr <-  unique(gsub("_","_chr",unique(unite(sseek, newcol, c(donor_unique_id,Chr), remove=FALSE,sep="_")$newcol)))
SReg_chr <-  unique(unite(SReg_results_in, newcol, c(duid,chrom), remove=FALSE,sep="_")$newcol)
SP_chr <- unique(unite(SP_chromo[which(SP_chromo$density < 5e-6),], newcol, c(sampleid,chrom), remove=FALSE,sep="_")$newcol)



# Fin the overlap between the three methods and plot a venn diagram

ll<- list(
    SReg_chr=SReg_chr,
    SP_chr=SP_chr,
    SSeek_chr=SSeek_chr)

venn.diagram(ll, filename ="~/tmp/shatterproof/venn_3methods_tmp.png",filetype="png")

hist(log10(sregs.bed$end-sregs.bed$start),xlim=c(4,9),breaks=100)
hist(log10(SSeek.bed$end-SSeek.bed$start),xlim=c(4,9),breaks=100)
hist(log10(SP_chromo$end-SP_chromo$start),xlim=c(4,9))

## evaluate the performace
# The total number of tests is defined by the number of samples in the gold standard times chromosomes (23):
a<- SReg_chr
b<- SP_chr

a<- SReg_chr
b<- SSeek_chr

a <- SP_chr
b<- SSeek_chr


tot <- length(unique(pcawg_chromo$donor_unique_id))*23
TP <- length(intersect(a,b))
FP <- length(setdiff(b,a))
FN <- length(setdiff(a,b))
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

### Evaluate the correlation of hot-spots

binsgr <- extract.bins(shreg_breast_pcawg)

sseek.mock <- sproof.mock <- sregs.mock <- list()
histoList <- names(which(sapply(unique(SSeek_results$histology_abbreviation),function(i) length(which(ali2hist == i))) > 20))

for(histology in histoList){
    
message(histology)
histoSamples <- ali2dui[names(which(ali2hist == histology))]

# create a HBD matrix from ShatterSeek bed file and save into a chromo.regs object

ss_hbd <- bed2chromo.reg(SSeek.bed[which(SSeek.bed$sampleid %in% histoSamples),],bingr)
sseek.mock[[histology]] <- new("chromo.regs",
                      high.density.regions=ss_hbd,
                      high.density.regions.hc=ss_hbd,
                      cnv=cnv,
                      svc=svc)
# create a HBD matrix from ShatterProof bed file and save into a chromo.regs object
sp_hbd <- bed2chromo.reg(SP_chromo[which(SP_chromo$sampleid %in% histoSamples),],bingr)
sproof.mock[[histology]] <- new("chromo.regs",
                   high.density.regions=sp_hbd,
                   high.density.regions.hc=sp_hbd,
                   cnv=cnv,
                   svc=svc)

# create a HBD matrix from ShatterProof bed file and save into a chromo.regs object
sr_hbd <- shreg_all_pcawg@high.density.regions.hc[intersect(rownames(shreg_all_pcawg@high.density.regions.hc),names(histoSamples)),]
sregs.mock[[histology]] <- new("chromo.regs",
                   high.density.regions=sr_hbd,
                   high.density.regions.hc=sr_hbd,
                   cnv=cnv,
                   svc=svc)
}


hh <- "Breast-AdenoCA"
# establish hot spot frequency cut off
par(mfrow=c(3,1))
shreg.fpt <- freq.p.test(sregs.mock[[hh]]@high.density.regions.hc)
sseek.fpt <- freq.p.test(sseek.mock[[hh]]@high.density.regions.hc)
sproof.fpt <- freq.p.test(sproof.mock[[hh]]@high.density.regions.hc)

# generate maps
par(mfrow=c(3,1))
shattered.map.plot(sregs.mock[[hh]],freq.cut = shreg.fpt@freq.cut)
shattered.map.plot(sseek.mock[[hh]], freq.cut = sseek.fpt@freq.cut)
shattered.map.plot(sproof.mock[[hh]], freq.cut = sproof.fpt@freq.cut)

# generate maps
hs_sregs <- hot.spot.samples(sregs.mock[[hh]],freq.cut = shreg.fpt@freq.cut)
hs_sseek <- hot.spot.samples(sseek.mock[[hh]],freq.cut = sseek.fpt@freq.cut)
hs_sproof <- hot.spot.samples(sproof.mock[[hh]],freq.cut = sproof.fpt@freq.cut)



shregfreq <- apply(sregs.mock[[hh]]@high.density.regions.hc,2,sum)
sseekfreq <- apply(sseek.mock[[hh]]@high.density.regions.hc,2,sum)
sprooffreq <- apply(sproof.mock[[hh]]@high.density.regions.hc,2,sum)


par(mfrow=c(1,3),mar=c(4,4,1,1))
ct1 <- cor.test(shregfreq,sseekfreq,method="pearson")
smoothScatter(shregfreq,sseekfreq,las=1)
points(shregfreq,sseekfreq,pch=19,cex=.2)
rect(shreg.fpt@freq.cut,sseek.fpt@freq.cut,100,100,border=NA,col=rgb(0.3, 0.3, 0.3, alpha = 0.3))
abline(v=shreg.fpt@freq.cut,lty=3,lwd=.5)
abline(h=sseek.fpt@freq.cut,lty=3,lwd=.5)
legend("topleft", paste("cor =",sprintf("%.2f",ct1$estimate),"\nP =",sprintf("%.2e",ct1$p.value)), bty='n')

ct2 <- cor.test(shregfreq,sprooffreq,method="pearson")
smoothScatter(shregfreq,sprooffreq,las=1)
points(shregfreq,sprooffreq,pch=19,cex=.2)
rect(shreg.fpt@freq.cut,sproof.fpt@freq.cut,100,100,border=NA,col=rgb(0.3, 0.3, 0.3, alpha = 0.3))
abline(v=shreg.fpt@freq.cut,lty=3,lwd=.5)
abline(h=sproof.fpt@freq.cut,lty=3,lwd=.5)
legend("topleft", paste("cor =",sprintf("%.2f",ct2$estimate),"\nP =",sprintf("%.2e",ct2$p.value)), bty='n')

ct3 <- cor.test(sseekfreq,sprooffreq,method="pearson")
smoothScatter(sseekfreq,sprooffreq,las=1)
points(sseekfreq,sprooffreq,pch=19,cex=.2)
rect(sseek.fpt@freq.cut,sproof.fpt@freq.cut,100,100,border=NA,col=rgb(0.3, 0.3, 0.3, alpha = 0.3))
abline(v=sseek.fpt@freq.cut,lty=3,lwd=.5)
abline(h=sproof.fpt@freq.cut,lty=3,lwd=.5)
legend("topleft", paste("cor =",sprintf("%.2f",ct3$estimate),"\nP =",sprintf("%.2e",ct3$p.value)), bty='n')

### functions:

extract.bins<- function(chromo.regs.obj){
         binNames <- colnames(chromo.regs.obj@high.density.regions)
         bindf <- remove.factors(data.frame(do.call(rbind,strsplit(binNames," "))))
         bindf[,2] <- as.numeric(bindf[,2])
         bindf[,3] <- as.numeric(bindf[,3])
         colnames(bindf) <- c("chrom","start","end")
         bingr <- with(bindf,GRanges(chrom, IRanges(start=start, end=end), name=binNames))
         return(bingr)
     }


bed2chromo.reg <- function(bed, bingr){

    inbedGR <- with(bed,GRanges(chrom, IRanges(start=start, end=end)))
    idNames <- unique(bed$sampleid)
    tot <- length(unique(bed$sampleid))*length(binNames)
    hbd <- matrix(rep(0,),ncol=length(binNames),nrow=length(idNames))
    colnames(hbd) <- bingr$name
    rownames(hbd) <- idNames
    
    hits <- GenomicAlignments::findOverlaps(bingr,inbedGR)
    
    for(i in idNames){
        hbdbins_i <- bingr$name[queryHits(hits)[which(bed$sampleid[subjectHits(hits)] == i)]]
        hbd[i,hbdbins_i] <- 1
    }
    return(hbd)
}

