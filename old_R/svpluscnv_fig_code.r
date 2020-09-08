library(svpluscnv)

setwd("<path-to-local-folder>/svpluscnv_doc_code/")

# data:
# CCLE BREAST data contains CNVs derived from SNP arrays and SVs derived from WGS
# download from DepMap (https://depmap.org/portal/download/)
load("data/brca_ccle.rda",verbose=TRUE)

# data:
# TCGA Breast-AdenoCA data contains CNVs and SVs derived from WGS
# download from PCAWG (https://dcc.icgc.org/releases/PCAWG)
load("data/brca_tcga.rda",verbose=TRUE)

# data:
# PCAWG all tumor types data contains CNVs and SVs derived from WGS
# download from GDC legacy portal (https://portal.gdc.cancer.gov/legacy-archive/
load("data/pcawg.rda",verbose=TRUE)
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

cnv.break.annot(cnv_brca_ccle,)
svc.break.annot(svc_brca_ccle)

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
sv.model.view(svc_brca_ccle, cnv_brca_ccle, chr, start, stop, sampleids=gene_samples, 
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

#### Figure Suppl S4: SHATTERED.TREGIONS

# run shattered.regions for each dataset
shreg_ccle <- shattered.regions(cnv_brca_ccle,svc_brca_ccle,fc.pct = 0.1, interleaved.cut = 0.33, verbose=FALSE)
shreg_tcga <- shattered.regions.cnv(cnv_brca_ccle, fc.pct = 0.1, verbose=FALSE)
shreg_pcawg <- shattered.regions(cnv_brca_pcawg,svc_brca_pcawg,fc.pct = 0.1, interleaved.cut = 0.33, verbose=FALSE)

# set a seed to fix randomization
set.seed=1234

# evaluate hotspots and plot shattered region maps for each dataset
pdf("figures/suppl_fig_S4.pdf",height=6,width=8)

layout(matrix(c(1,2,3,4,5,6),3,2,
              byrow = TRUE),widths = c(2,6))
par(mar=c(4,4,3,1))
fpt_ccle <- freq.p.test(shreg_ccle@high.density.regions.hc)
par(mar=c(3,4,3,4))
shattered.map.plot(shreg_ccle,freq.cut = fpt_ccle@freq.cut)
text(4e8,fpt_ccle@freq.cut+0.5,"fdr < 0.05",cex=1.1)

par(mar=c(4,4,3,1))
fpt_tcga <- freq.p.test(shreg_tcga@high.density.regions.hc)
par(mar=c(3,4,3,4))
shattered.map.plot(shreg_tcga,freq.cut = fpt_tcga@freq.cut)
text(4e8,fpt_tcga@freq.cut+0.6,"fdr < 0.05",cex=1.1)

par(mar=c(4,4,3,1))
fpt_pcawg <- freq.p.test(shreg_pcawg@high.density.regions.hc)
par(mar=c(3,4,3,4))
shattered.map.plot(shreg_pcawg,freq.cut = fpt_pcawg@freq.cut)
text(4e8,fpt_pcawg@freq.cut+1,"fdr < 0.05",cex=1.1)

dev.off()


### FIGURE 1

# layout.show(4)
pdf("figures/Figure1.pdf",height=8,width=8)
layout(matrix(c(1,1,2,3,4,4),3,2,byrow = TRUE),heights = c(1.5,3,1.5))
#1
par(mar=c(1,4,1,4))
cnv.freq(cnv_brca_ccle,fc.pct = 0.3, ploidy = TRUE)
#2
par(mar=c(1,1,4,1),font.sub=3)
circ.chromo.plot(shreg_ccle,sample.id = "MDAMB134VI_BREAST")
legend("topleft",c("DUP","DEL","INV","TRA"),lwd=1,
       col=c("red","blue","gold","black"),
       title= expression(bold(underline("SVC"))),bty='n')
legend("topright",c("Gain","2n","Loss"),lwd=1,
       col=c("red","black","blue"),
       title=expression(bold(underline("CNV"))),bty='n')
#3
par(mar=c(4,4,4,1))
fpt_ccle <- freq.p.test(shreg_ccle@high.density.regions.hc)
#4
par(mar=c(3,4,1,4))
shattered.map.plot(shreg_ccle,freq.cut = fpt_ccle@freq.cut)
text(4e8,fpt_ccle@freq.cut+1,"fdr < 0.05",cex=1.1)
dev.off()



