###
# Author: Gonzalo Lopez, PhD
# e-mail: gonzolgarcia@gmail.com; gonzalo.lopezgarcia@mssm.edu

### prepare PCAWG dataset to run ShatterProof algorithm 
# (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-78)


require(data.table)
require(tidyr)
require(GenomicRanges)
require(taRifx)

# move to local dir
setwd("<local-path-to>/svpluscnv_doc_code/")

# data:
# PCAWG all tumor types data contains CNVs and SVs derived from WGS
# download from GDC legacy portal (https://portal.gdc.cancer.gov/legacy-archive/
load("data/pcawg.rda",verbose=TRUE)

# prepare input for ShatterProof; only CN and TRA filetypes available

smaple_id_list <- intersect(cnv_pcawg@data$sample,svc_pcawg@data$sample)
write.table(smaple_id_list,file="shatterproof/list_of_samples.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

dir.create("shatterproof/cn")
options(scipen=999)
for(i in smaple_id_list){ 
    message(i)
    segment_cn_i <- segment_cn[which(segment_cn$sample == i)][,c(2,3,4,6,5)]
    colnames(segment_cn_i) <- c("#chr","start","end","number","quality")
    segment_cn_i$quality[] <- "."
    segment_cn_i$`#chr` <- gsub("chr","", segment_cn_i$`#chr`)
    
    write.table(segment_cn_i,file=paste("shatterproof/cn/",i,".spc",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
}

tra_class <- which(svc_pcawg@data$svclass == "TRA")
non_tra_class <- which(svc_pcawg@data$svclass != "TRA")
mb2plus <- intersect(non_tra_class,which(svc_pcawg@data$pos2 -svc_pcawg@data$pos1 > 2e6) )
svc_pcawg_TRA <- svc_pcawg@data[sort(c(tra_class,mb2plus))]
dir.create("shatterproof/tra")
for(i in smaple_id_list){ 
    svc_pcawg_TRA_i <- svc_pcawg_TRA[which(svc_pcawg_TRA$sample == i)][,c(2,3,3,5,6,6,7)]
    colnames(svc_pcawg_TRA_i) <- c("#chr1","start","end","chr2","start","end","quality")
    svc_pcawg_TRA_i$quality[] <- "."
    svc_pcawg_TRA_i$`#chr1` <- gsub("chr","", svc_pcawg_TRA_i$`#chr1`)
    svc_pcawg_TRA_i$`chr2` <- gsub("chr","", svc_pcawg_TRA_i$`chr2`)
    write.table(svc_pcawg_TRA_i,file=paste("shatterproof/tra/",i,".spt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    message(paste(i,nrow(svc_pcawg_TRA_i)))
}
options(scipen=0)

dir.create("shatterproof/out")

# Run shatterproof using perk scrips shatterproof/run_shatterproof.pl:
# 
# data:

sp_filedata <- list()
sp_chromo_files <- dir("shatterproof/out//")
for(i in sp_chromo_files){
    SPfileName <- paste("shatterproof/out/",i,"/suspect_regions/suspect_regions.yml",sep="")
    if(file.exists(SPfileName)){
        nLines <- as.numeric(strsplit(system(paste("wc",SPfileName), intern = TRUE)," +")[[1]][2])
        message(paste(i,nLines,sep=": "))
        dat <-  remove.factors(read.delim(SPfileName,sep=";",header=FALSE))
        regions <- grep("chromosome:|start:|end:",dat$V1, perl=T,value=T)
        df <- data.frame(t(matrix(regions,nrow = 4,ncol = length(regions)/4)))
        df[,1] <- gsub("chromosome:\t","chr",df[,1])
        df[,2] <- gsub("start:\t\t","",df[,2])
        df[,3] <- gsub("end:\t\t","",df[,3])
        df[,4] <- gsub("mutation_density_of_chromosome:\t","",df[,4])
        for(j in 2:4) df[,j] <- as.numeric(df[,j])
        colnames(df) <- c("chrom","start","end","density")
        sampleid <- rep(ali2dui[i],nrow(df))
        name2 <- rep(i,nrow(df))
        sp_filedata[[i]] <-data.frame(sampleid,df,name2)
    }
}
SP_chromo <- remove.factors(do.call(rbind,sp_filedata))

## save all likely chromothriptic regions for further evaluation
save(SP_chromo, file="data/shatterproof_pcawg_chromo_results.rda")



