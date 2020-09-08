library(svpluscnv)
library(tidyr)
library(data.table)
library(VennDiagram)

setwd("<path-to-local-folder>/svpluscnv_doc_code/")
setwd("~/Box Sync/git/svpluscnv_doc_code/")

# data:
# CCLE BREAST data contains CNVs derived from SNP arrays and SVs derived from WGS
# download from DepMap (https://depmap.org/portal/download/)
load("data/brca_ccle.rda",verbose=TRUE)

# data:
# TCGA Breast-AdenoCA data contains CNVs from the GDC legacy archive
# download from GDC legacy portal (https://portal.gdc.cancer.gov/legacy-archive/
load("data/brca_tcga.rda",verbose=TRUE)

# data:
# PCAWG all tumor types data contains CNVs and SVs derived from WGS
# download from PCAWG (https://dcc.icgc.org/releases/PCAWG)

load("data/pcawg.rda",verbose=TRUE)
# extract Breast Cancer data from PCAWG
brca_samples_pcawg <- names(which(ali2hist == "Breast-AdenoCA"))
cnv_brca_pcawg <- cnv_pcawg; cnv_brca_pcawg@data <- cnv_brca_pcawg@data[which(cnv_brca_pcawg@data$sample %in% brca_samples_pcawg)]
svc_brca_pcawg <- svc_pcawg; svc_brca_pcawg@data <- svc_brca_pcawg@data[which(svc_brca_pcawg@data$sample %in% brca_samples_pcawg)]


# Figure Suppl 1
pdf("figures/suppl_fig_S1.pdf",height=5,width=6)
par(mfrow=c(3,1),mar=c(3,4,1,4))
cnv.freq(cnv_brca_ccle,fc.pct = 0.3, ploidy = TRUE)
cnv.freq(cnv_brca_tcga,fc.pct = 0.3, ploidy = TRUE)
cnv.freq(cnv_brca_pcawg,fc.pct = 0.3, ploidy = TRUE)
dev.off()

#Figure Suppl 2

# load a list of cancer genes (cosmic cancer census release v90)
load("data/cosmic_census_v90.rda",verbose=TRUE)

## Identification of breakpoints overlapping known genes; for CCLE we use bothe CNV and SV and intersect the results
cnv_annot_ccle <- cnv.break.annot(cnv_brca_ccle, 
                               fc.pct = 0, 
                               genome.v="hg19",
                               clean.brk = 8)
svc_annot_ccle <- svc.break.annot(svc_brca_ccle, 
                               svc.seg.size = 200000, 
                               genome.v="hg19")
# combine CNV and SVC and filter by cosmic genes
disruptSamplesCCLE <- merge2lists(cnv_annot_ccle@disruptSamples,svc_annot_ccle@disruptSamples, fun="intersect")
disruptSamplesCCLE_cosmic <- disruptSamplesCCLE[intersect(names(disruptSamplesCCLE),cosmic_cancer_census_v90$Gene.Symbol)]
ccle_plot_genes <- rev(sort(unlist(lapply(disruptSamplesCCLE_cosmic,length)),decreasing=T)[1:20])

# Identification of breakpoints overlapping known genes; for TCGA dataset we only use CNV 

cnv_annot_tcga <- cnv.break.annot(cnv_brca_tcga, 
                                  fc.pct = 0, 
                                  genome.v="hg19", 
                                  clean.brk = 8)
# filter by cosmic genes
disruptSamplesTCGA <- cnv_annot_tcga@disruptSamples
disruptSamplesTCGA_cosmic <- disruptSamplesTCGA[intersect(names(disruptSamplesTCGA),cosmic_cancer_census_v90$Gene.Symbol)]
tcga_plot_genes <- rev(sort(unlist(lapply(disruptSamplesTCGA_cosmic,length)),decreasing=T)[1:20])


## Identification of breakpoints overlapping known genes; for PCAWG we use bothe CNV and SV and intersect the results 
cnv_annot_pcawg <- cnv.break.annot(cnv_brca_pcawg, 
                                   fc.pct = 0,
                                   genome.v="hg19", 
                                   clean.brk = 8)
svc_annot_pcawg <- svc.break.annot(svc_brca_pcawg, 
                                   svc.seg.size = 200000, 
                                   genome.v="hg19", 
                                   verbose=FALSE)
# combine CNV and SVC and filter by cosmic genes
disruptSamplesPCAWG <- merge2lists(cnv_annot_pcawg@disruptSamples, svc_annot_pcawg@disruptSamples, fun="intersect")
disruptSamplesPCAWG_cosmic <- disruptSamplesPCAWG[intersect(names(disruptSamplesPCAWG),cosmic_cancer_census_v90$Gene.Symbol)]
pcawg_plot_genes <- rev(sort(unlist(lapply(disruptSamplesPCAWG_cosmic,length)),decreasing=T)[1:20])


## obtain gene coordinates for ploting 
gene <- "FHIT"
df <- gene.track.view(symbol = gene, plot=FALSE, genome.v = "hg19")@data
start <- min(df$txStart) - 50000;  stop <- max(df$txEnd) + 300000;  chr <- df$chrom[1]
gene_samples <- sort(disruptSamplesCCLE[[gene]])

# Obtain genes in the top 20 in two datasets
repl2 <- names(which(sort(table(c(names(ccle_plot_genes),names(tcga_plot_genes),names(pcawg_plot_genes)))) == 2))
# Obtain genes in the top 20 in three datasets
repl3 <- names(which(sort(table(c(names(ccle_plot_genes),names(tcga_plot_genes),names(pcawg_plot_genes)))) == 3))
# Create a gene color vector for each dataset based on replication across datasets
ccle_color <- tcga_color <- pcawg_color <- rep("grey",20)
names(ccle_color) <- names(ccle_plot_genes)
names(tcga_color) <- names(tcga_plot_genes)
names(pcawg_color) <- names(pcawg_plot_genes)
ccle_color[repl2] <- "lightblue"; ccle_color[repl3] <- "blue"; 
tcga_color[repl2] <- "lightblue"; tcga_color[repl3] <- "blue"
pcawg_color[repl2] <- "lightblue"; pcawg_color[repl3] <- "blue"

# plot figure suppl S2
pdf("figures/suppl_fig_S2.pdf",height=8,width=6)
layout(matrix(c(1,2,3,4,4,4,5,5,5),3,3,
              byrow = TRUE),heights = c(5.5,4.5,3.5))

par(mar=c(4,5,1,1))
barplot(ccle_plot_genes,
        horiz=T,las=1,cex.names = 1, xlab="#CCLE samples",
        col=ccle_color,border="NA")
barplot(tcga_plot_genes,
        horiz=T,las=1,cex.names = 1, xlab="#TCGA samples",
        col=tcga_color,border="NA")
barplot(pcawg_plot_genes,
        horiz=T,las=1,cex.names = 1, xlab="#PCAWG samples",
        col=pcawg_color,border="NA")

par(mar=c(1,10,0,1))
sv.model.view(cnv_brca_ccle, svc_brca_ccle, chr, start, stop, sampleids=gene_samples, 
              addlegend = 'both', addtext=c("TRA"), cnvlim = c(-2,2), 
              cex=.7,cex.text =.8, summary = FALSE,cex.legend = 0.7)
par(mar=c(3,10,1,1))
gene.track.view(chr=chr ,start=start, stop=stop, addtext=TRUE, cex.text=1, 
                summary = FALSE,cex=1.5)

dev.off()


#### Figure SUPPL S3

# collect breakpoints for CNV and SVC
cnv_breaks_ccle <- cnv.breaks(cnv_brca_ccle)
svc_breaks_ccle <- svc.breaks(svc_brca_ccle)
# match both breakpoint sources and plot
pdf("figures/suppl_fig_S3.pdf",height=5,width=6)
par(mar=c(5,5,2,2))
match.breaks(cnv_breaks_ccle,svc_breaks_ccle, maxgap = 10000)
dev.off()

#### #### #### #### #### #### #### #### #### #### #### #### 
#### Figure Suppl S4: SHATTERED.TREGIONS
# set a seed to fix randomization
set.seed=1234

# run shattered.regions for each dataset
## We defiune chromosome list and coordinates in order to be able to compare across different runs
chrlist = paste("chr",c(1:22,"X"),sep="")
chr.lim.hg19 <- d3gb.chr.lim(genome.v = "hg19")
#chr.lim.hg19 <- chr.lim.hg19[which(chr.lim.hg19$chrom %in% chrlist)]
chr.lim.hg38 <- d3gb.chr.lim(genome.v = "hg38")


shreg_ccle <- shattered.regions(cnv_brca_ccle,svc_brca_ccle, 
                                fc.pct = 0.1, interleaved.cut = 0.33, 
                                chrlist=chrlist, chr.lim = chr.lim.hg19,
                                verbose=FALSE)
shreg_tcga <- shattered.regions.cnv(cnv_brca_tcga, 
                                fc.pct = 0.1,
                                chrlist = chrlist, chr.lim = chr.lim.hg19,
                                verbose=FALSE)
shreg_pcawg <- shattered.regions(cnv_brca_pcawg,svc_brca_pcawg,
                                fc.pct = 0.1, interleaved.cut = 0.33, 
                                chrlist = chrlist, chr.lim = chr.lim.hg38,
                                verbose=FALSE)


### NEW FIGURE 1

## CHROMOTHRIPSIS analysis of chromothripsis based on ShatterProof, ShatterSeek and shattered.regions (svpluscnv);
# comparison across methods

## chromothripsis data from PCAWG
# Chromothripsis survey from PCAWG table contains info about complex chromosomal rearrangements in PCAWG samples
# download from manuscript (https://doi.org/10.1038/s41588-019-0576-7) based on ShatterSeek algorithm (https://github.com/parklab/ShatterSeek)
# suppl table 1 (https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0576-7/MediaObjects/41588_2019_576_MOESM3_ESM.xlsx)

load("data/pcawg_chromo.rda",verbose=TRUE)
# reformat to bed format
SSeek_results <- pcawg_chromo[which(pcawg_chromo$chromo_label %in% c("High confidence","Low confidence") ),]
SSeek_results_breast <- SSeek_results[which(SSeek_results$histology_abbreviation == "Breast-AdenoCA" ),]
SSeek.bed <- as.data.frame(SSeek_results[,c("donor_unique_id","Chr","Start","End")])
SSeek.bed[,"Chr"] <- paste("chr",SSeek.bed[,"Chr"] ,sep="")
colnames(SSeek.bed) <- c("sampleid","chrom","start","end")

# load data from 
# Generated using code available at ~/shatterproof/prepare_shatterproof.r and ~/shatterproof/run_shatterproof.pl
load("data/shatterproof_pcawg_chromo_results.rda",verbose=TRUE)
SProof.bed <- SP_chromo[SP_chromo$sampleid %in% unique(pcawg_chromo$donor_unique_id),]


# run shattered.regions for all PCAWG samples
#set.seed(1234)
#shreg_all_pcawg <- shattered.regions(cnv_pcawg, svc_pcawg, fc.pct = 0.1, interleaved.cut = 0.33, verbose=TRUE)
#save(shreg_all_pcawg, file="data/shreg_all_pcawg.rda")
load("data/shreg_all_pcawg.rda",verbose=TRUE)

pdf("figures/Figure1_A_CircosPlot.pdf",width=5,height=5)
circ.chromo.plot(shreg_all_pcawg,sample.id = "fc7f8eeb-9c40-123e-e040-11ac0c484061",
                 high.conf = TRUE,print.name = FALSE)
dev.off()

# in order to link the sampleids we are going to use donor_unique_id field
spcid <- unname(unlist(sapply(names(shreg_all_pcawg@regions.summary), 
                              function(i) rep(i, nrow(shreg_all_pcawg@regions.summary[[i]])))))
duid <- ali2dui[spcid]  
SReg_results <- data.table(spcid,duid,do.call(rbind,shreg_all_pcawg@regions.summary))

# to make data comparable we restrict the analysis to samples included by the PCAWG study only (see https://doi.org/10.1038/s41588-019-0576-7 for filtering details) 
SReg_results_in <- SReg_results[which(SReg_results$duid %in% unique(pcawg_chromo$donor_unique_id)),]
SRegs.bed <- SReg_results_in[,2:5]   
colnames(SRegs.bed) <- c("sampleid","chrom","start","end")

# paste donor_id and chromosome from shattered regions from shattered.regions and PCAWG study
SSeek_chr <-  unique(unite(SSeek.bed, newcol, c(sampleid,chrom), remove=FALSE,sep="_")$newcol)
SRegs_chr <-  unique(unite(SRegs.bed, newcol, c(sampleid,chrom), remove=FALSE,sep="_")$newcol)
SProof_chr <- unique(unite(SProof.bed, newcol, c(sampleid,chrom), remove=FALSE,sep="_")$newcol)

# obtain overlap of shattered chromosomes across three methods
ll<- list(
        shattered.regions=SRegs_chr,
        ShatterProof=SProof_chr,
        ShatterSeek=SSeek_chr
)
venn.diagram(ll, filename ="figures/figure1A.png",filetype="png")

pdf("figures/suppementary_S5ABC.pdf",width=8,height = 3)
par(mfrow=c(1,3),mar=c(4,5,4,1))
hist(log10(SRegs.bed$end-SRegs.bed$start),xlim=c(4,9),breaks=50,xlab="log10(bp)",las=1,main="shattered.regions")
hist(log10(SSeek.bed$end-SSeek.bed$start),xlim=c(4,9),breaks=50,xlab="log10(bp)",las=1,main="ShatterSeek")
hist(log10(SProof.bed$end-SProof.bed$start),xlim=c(4,9),breaks=10,xlab="log10(bp)",las=1,main="ShatterProof")
dev.off()

#The total number of tests is defined by the number of samples in the gold standard times chromosomes (23):

accprerec <- function(a,b,tot){
        TP <- length(intersect(a,b))
        FP <- length(setdiff(a,b))
        FN <- length(setdiff(b,a))
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
        res <- c(Accuracy,Precision,Recall)
        names(res) <- c("Accuracy","Precision","Recall")
        return(res)
}

tot <- length(unique(pcawg_chromo$donor_unique_id))*23

ll<-list()
a<- SRegs_chr
b<- unique(c(SProof_chr,SSeek_chr))
ll[["sregs"]] <- accprerec(a,b,tot)

a<- SProof_chr
b<- unique(c(SSeek_chr,SRegs_chr))
ll[["sproof"]] <- accprerec(a,b,tot)

a <- SSeek_chr
b <- unique(c(SProof_chr,SRegs_chr))
ll[["sseek"]] <- accprerec(a,b,tot)

pdf("figures/figure1B.pdf",width=4.5,heigh=3.5)
barplot(t(do.call(rbind,ll)[,2:3]),beside = TRUE,las=1,names=c("shattered.reg","ShatterProof","ShatterSeek"),ylim=c(0,0.99))
grid(nx=NA, ny =NULL)
legend("topleft",c("Precission","Recall"),fill=c("grey30","grey80"),ncol=2,bty='n')
dev.off()


ll<-list()
a<- unique(SRegs_chr)
b<- unique(c(SSeek_chr))
ll[["sregs"]] <- accprerec(a,b,tot)
ab <- length(intersect(a,b))
x <- fisher.test(matrix(c(ab, length(a) -ab, length(b)-ab, tot -length(a) -length(b) + ab),2,2),alternative = "greater")
x$p.value

a<- unique(SProof_chr)
b<- unique(c(SSeek_chr))
ll[["sproof"]] <- accprerec(a,b,tot)
x <- fisher.test(matrix(c(ab, length(a) -ab, length(b)-ab, tot -length(a) -length(b) + ab),2,2),alternative = "greater")
x$p.value

a<- unique(SProof_chr)
b<- unique(SRegs_chr)
x <- fisher.test(matrix(c(ab, length(a) -ab, length(b)-ab, tot -length(a) -length(b) + ab),2,2),alternative = "greater")
x$p.value


pdf("figures/figure1Bbis.pdf",width=4.5,heigh=3.5)
barplot(t(do.call(rbind,ll)),beside = TRUE,las=1,names=c("shattered.reg","ShatterProof"),ylim=c(0,1))
grid(nx=NA, ny =NULL)
#legend("top",c("Accuracy","Precission","Recall"),fill=c("grey30","grey50","grey80"),ncol=3,bty='n')
dev.off()
pdf("figures/figure1B_barplot_legend.pdf",width=4.5,heigh=3.5)
plot(1,1,col="white",xaxt=',',yaxt='n',bty='n')
legend("top",c("Accuracy","Precission","Recall"),fill=c("grey30","grey50","grey80"),ncol=1,bty='n')
dev.off()

### Evaluate the correlation of hot-spots

binsgr <- extract.bins(shreg_all_pcawg)

sseek.mock <- sproof.mock <- sregs.mock <- list()
histoList <- names(which(sapply(unique(SSeek_results$histology_abbreviation),function(i) length(which(ali2hist == i))) > 20))

for(histology in histoList){
        
        message(histology)
        histoSamples <- ali2dui[names(which(ali2hist == histology))]
        
        # create a HBD matrix from ShatterSeek bed file and save into a chromo.regs object
        
        ss_hbd <- bed2chromo.reg(SSeek.bed[which(SSeek.bed$sampleid %in% histoSamples),],binsgr)
        sseek.mock[[histology]] <- new("chromo.regs",
                                       high.density.regions=ss_hbd,
                                       high.density.regions.hc=ss_hbd)
        # create a HBD matrix from ShatterProof bed file and save into a chromo.regs object
        sp_hbd <- bed2chromo.reg(SP_chromo[which(SP_chromo$sampleid %in% histoSamples),],binsgr)
        sproof.mock[[histology]] <- new("chromo.regs",
                                        high.density.regions=sp_hbd,
                                        high.density.regions.hc=sp_hbd)
        
        # create a HBD matrix from ShatterProof bed file and save into a chromo.regs object
        sr_hbd <- shreg_all_pcawg@high.density.regions.hc[intersect(rownames(shreg_all_pcawg@high.density.regions.hc),names(histoSamples)),]
        sregs.mock[[histology]] <- new("chromo.regs",
                                       high.density.regions=sr_hbd,
                                       high.density.regions.hc=sr_hbd)
}

ss_hbd <- bed2chromo.reg(SSeek.bed[which(SSeek.bed$sampleid %in% histoSamples),],binsgr)

cl_hbd <- shreg_ccle@high.density.regions.hc
sregs.ccle <- new("chromo.regs",
                               high.density.regions=cl_hbd,
                               high.density.regions.hc=cl_hbd)

colnames(shreg_ccle@high.density.regions)[1:10]
colnames(shreg_all_pcawg@high.density.regions.hc)[1:10]

### Shattered maps and hot-spot






hh <- "Breast-AdenoCA"
# establish hot spot frequency cut off
pdf("figures/suppl_fig_S4.pdf",height=5,width=7.5)
par(mfrow=c(2,3),mar=c(5,5,3,1))
shreg.fpt <- freq.p.test(sregs.mock[[hh]]@high.density.regions.hc)
title("PCAWG: shattered.resions",cex.main=1)
sseek.fpt <- freq.p.test(sseek.mock[[hh]]@high.density.regions.hc)
title("PCAWG: ShatterSeek",cex.main=1)
sproof.fpt <- freq.p.test(sproof.mock[[hh]]@high.density.regions.hc)
title("PCAWG: ShatterProof",cex.main=1)
fpt_ccle <- freq.p.test(shreg_ccle@high.density.regions.hc)
title("CCLE: shattered.resions",cex.main=1)
fpt_tcga <- freq.p.test(shreg_tcga@high.density.regions.hc)
title("TCGA: shattered.resions",cex.main=1)
dev.off()

# plot only chromosome 1-22 and X

pdf("figures/Figure1_sregs.pdf",width=8,height=0.7)
par(mar=c(0,4,0,4))
shattered.map.plot(sregs.mock[[hh]],freq.cut = shreg.fpt@freq.cut,add.legend = FALSE, chrlist=chrlist)
dev.off()

pdf("figures/Figure1_sseek.pdf",width=8,height=0.7)
par(mar=c(0,4,0,4))
shattered.map.plot(sseek.mock[[hh]], freq.cut = sseek.fpt@freq.cut,add.legend = FALSE, chrlist=chrlist)
dev.off()

pdf("figures/Figure1_sproof.pdf",width=8,height=0.7)
par(mar=c(0,4,0,4))
shattered.map.plot(sproof.mock[[hh]], freq.cut = sproof.fpt@freq.cut,add.legend = FALSE, chrlist=chrlist)
dev.off()

pdf("figures/Figure1_ccle.pdf",width=8,height=0.7)
par(mar=c(0,4,0,4))
shattered.map.plot(shreg_ccle,freq.cut = fpt_ccle@freq.cut,add.legend = FALSE, chrlist=chrlist)
dev.off()

pdf("figures/Figure1_tcga.pdf",width=8,height=0.7)
par(mar=c(0,4,0,4))
shattered.map.plot(shreg_tcga,freq.cut = fpt_tcga@freq.cut, add.legend = FALSE, chrlist=chrlist)
dev.off()

# generate maps
hs_sregs <- hot.spot.samples(sregs.mock[[hh]],freq.cut = shreg.fpt@freq.cut)
hs_sseek <- hot.spot.samples(sseek.mock[[hh]],freq.cut = sseek.fpt@freq.cut)
hs_sproof <- hot.spot.samples(sproof.mock[[hh]],freq.cut = sproof.fpt@freq.cut)
hs_sccle <- hot.spot.samples(shreg_ccle,freq.cut = fpt_ccle@freq.cut)
hs_stcga <- hot.spot.samples(shreg_tcga,freq.cut = fpt_tcga@freq.cut)



shregfreq <- apply(sregs.mock[[hh]]@high.density.regions.hc,2,sum)
sseekfreq <- apply(sseek.mock[[hh]]@high.density.regions.hc,2,sum)
sprooffreq <- apply(sproof.mock[[hh]]@high.density.regions.hc,2,sum)
scclefreq <- apply(shreg_ccle@high.density.regions.hc,2,sum)
stcgafreq <- apply(shreg_tcga@high.density.regions.hc,2,sum)

mycol <- rgb(255, 255, 255, max = 255, alpha = 155)

pdf("Figures/corrplots_v3.pdf",width=6,height=6)
par(mfrow=c(2,2),mar=c(4,4,1,1))
ct1 <- cor.test(shregfreq,sseekfreq,method="pearson")
smoothScatter(shregfreq,sseekfreq,las=1)
points(shregfreq,sseekfreq,pch=19,cex=.3)
rect(shreg.fpt@freq.cut,sseek.fpt@freq.cut,100,100,border=NA,col=rgb(0.3, 0.3, 0.3, alpha = 0.3))
abline(v=shreg.fpt@freq.cut,lty=3,lwd=.3)
abline(h=sseek.fpt@freq.cut,lty=3,lwd=.3)
legend("bottomright", paste("cor =",sprintf("%.2f",ct1$estimate),"\nP =",sprintf("%.2e",ct1$p.value)), bty='n',cex=1.1)

ct2 <- cor.test(shregfreq,sprooffreq,method="pearson")
smoothScatter(shregfreq,sprooffreq,las=1)
points(shregfreq,sprooffreq,pch=19,cex=.3)
rect(shreg.fpt@freq.cut,sproof.fpt@freq.cut,100,100,border=NA,col=rgb(0.3, 0.3, 0.3, alpha = 0.3))
abline(v=shreg.fpt@freq.cut,lty=3,lwd=.3)
abline(h=sproof.fpt@freq.cut,lty=3,lwd=.3)
legend("bottomright", paste("cor =",sprintf("%.2f",ct2$estimate),"\nP =",sprintf("%.2e",ct2$p.value)), bty='n',cex=1.1)

ct3 <- cor.test(sseekfreq,sprooffreq,method="pearson")
smoothScatter(sseekfreq,sprooffreq,las=1)
points(sseekfreq,sprooffreq,pch=19,cex=.2)
rect(sseek.fpt@freq.cut,sproof.fpt@freq.cut,100,100,border=NA,col=rgb(0.3, 0.3, 0.3, alpha = 0.3))
abline(v=sseek.fpt@freq.cut,lty=3,lwd=.5)
abline(h=sproof.fpt@freq.cut,lty=3,lwd=.5)
legend("bottomright", paste("cor =",sprintf("%.2f",ct3$estimate),"\nP =",sprintf("%.2e",ct3$p.value)), bty='n',cex=1.1)

ct4 <- cor.test(scclefreq,stcgafreq,method="pearson")
smoothScatter(scclefreq,stcgafreq,las=1,cex.lab=1.5)
points(scclefreq,stcgafreq,pch=19,cex=.3)
rect(fpt_ccle@freq.cut,fpt_tcga@freq.cut,100,100,border=NA,col=rgb(0.3, 0.3, 0.3, alpha = 0.3))
abline(v=fpt_ccle@freq.cut,lty=3,lwd=.3)
abline(h=fpt_tcga@freq.cut,lty=3,lwd=.3)
legend("bottomright", paste("cor =",sprintf("%.2f",ct4$estimate),"\nP =",sprintf("%.2e",ct4$p.value)), box.col=mycol,cex=1.1,bg=mycol)


dev.off()




breastchromo <- names(which(ali2hist == "Breast-AdenoCA"))[which(names(which(ali2hist == "Breast-AdenoCA")) %in% names(shreg_all_pcawg@regions.summary))]

shreg_all_pcawg@regions.summary[breastchromo]

#sample.id <- "fc8130df-3147-3e94-e040-11ac0d485df8"
#sample.id <- "f393bb05-53c2-f80a-e040-11ac0d484528"
#sample.id <- "f393bafd-1baa-e5f4-e040-11ac0d48450b"
#sample.id <- "fc8130e0-0b9c-bbc9-e040-11ac0c483266"
#sample.id <- "f7fdda4f-7bf7-ede7-e040-11ac0c486e57"
#sample.id <- "74039acd-5aca-4c65-818c-3b577d295be0"


sample.id <- "fc7f8eeb-9c40-123e-e040-11ac0c484061"
pdf("figures/Figure1A_R1.pdf",width=8,height=8)
circ.chromo.plot(shreg_all_pcawg, sample.id = sample.id, high.conf = TRUE,print.name=FALSE)
dev.off()



