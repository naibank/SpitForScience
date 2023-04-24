### the first line only need to run once to install the GSBurden package
#install.packages(c("rstudioapi", "GenomicRanges", "BiasedUrn", "MASS", "ordinal"))

library(GenomicRanges)
library(GSBurden)
library(ggplot2)

#################################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sampleinfo <- readxl::read_excel("../phenotypes/Supplementary Table S1.v13.xlsx", sheet = 1, 
                                 col_types = c(rep("text", 11), rep("numeric", 29)))
sampleinfo <- sampleinfo[sampleinfo$ethnicity == "EAS", ]
sampleinfo$Plate <- paste0("Plate",gsub("Plate\\ ", "", sampleinfo$Plate))
sampleinfo$Plate <- gsub("Plate\\ ", "", sampleinfo$Plate)
gene.in <- read.delim("../../../ReferenceData/geneInfo2019/hg19_refFlat_exon.tsv", stringsAsFactors = F)
gene.in <- gene.in[which(gene.in$type_of_gene == "protein-coding"), ]
genes <- GeneAnnotation(gene.in$enzid, gene.in$chr, gene.in$start, gene.in$end, gene.in$genesymbol)
pli <- read.delim("../../../ReferenceData/constraint/gnomad.v2.1.1.lof_metrics.by_gene_March2020_hg37.txt", stringsAsFactors = F)
# genes <- genes[genes$gsymbol %in% pli$gene[which(pli$pLI > 0.9)], ]

### read cnv data
cnv.osc1 <- read.delim("../cnvs/Fragmentation_Merge.cnv.annotation20190807.Annot.txt", stringsAsFactors = F)
cnv.osc1 <- cnv.osc1[cnv.osc1$stringentFlag == "stringent", ]
cnv.osc1 <- subset(cnv.osc1, size >= 10000 & probeNum >= 5 &
                     unclean < 70 & segdup < 70 &  repeats < 70 &
                     DGV.OnePrcnt.Clusters == "No.Cluster(s)=0") 
cnv.osc1$Against.East.Asian.Int.Freq <- cnv.osc1$Against.East.Asian.Int.Freq - min(cnv.osc1$Against.East.Asian.Int.Freq, na.rm = T)
cnv.osc1$Against.Eur.Int.Freq <- cnv.osc1$Against.Eur.Int.Freq - min(cnv.osc1$Against.Eur.Int.Freq, na.rm = T)

cnv.osc1$max_freq <- pmax(cnv.osc1$Against.All.TCAG.Conts.Freq, pmax(cnv.osc1$Against.Eur.TCAG.Cont.Freq, 
                                                                     pmax(cnv.osc1$Against.Chinese.Conts.Freq,
                                                                          pmax(cnv.osc1$Against.Eur.Int.Freq, cnv.osc1$Against.East.Asian.Int.Freq, na.rm = T), na.rm = T), na.rm = T), na.rm = T)

cnv.osc2 <- read.delim("../cnvs/sfs3141passQC.TAGID_gender.cnv_annotation20190807.v2.txt", stringsAsFactors = F)
cnv.osc2 <- cnv.osc2[cnv.osc2$stringentFlag == "stringent", ]
cnv.osc2 <- subset(cnv.osc2, size >= 10000 & probeNum >= 5 &
                     unclean < 70 & segdup < 70 &  repeats < 70 &
                     DGV.OnePrcnt.Clusters == "No.Cluster(s)=0")
cnv.osc2$Internal.Freq <- cnv.osc2$Internal.Freq-min(cnv.osc2$Internal.Freq, na.rm = T)
cnv.osc2$max_freq <- pmax(cnv.osc2$Against.TCAG.Conts.All.Freq, cnv.osc2$Internal.Freq, na.rm = T)

feat <- intersect(names(cnv.osc1), names(cnv.osc2))
cnv <- rbind(cnv.osc1[, feat], cnv.osc2[, feat])
cnv$max_freq[is.na(cnv$max_freq)] <- 0

samples_with_CNV3mb <- cnv$SID[cnv$size >= 3000000]
writeLines(samples_with_CNV3mb, "samples_with_CNV3mb.txt")

### read snv data
outliers <- readLines("../qcscripts/outliers.inCNVs.txt") #, readLines("../qcscripts/outliers.inAll.LoF.txt"))
# snv <- data.table::fread("../snvs/OSC1_OSC2.All.SNVs.hg19.Eur.EastAsian.10prct.rm_too_high_internal.txt", data.table = F)
# snv <- snv[which((snv$LOF == "LOF" & snv$gnomAD_oe_lof_upper < 0.35) |
#                    (snv$MSN == "DMSN" & snv$LOF != "LOF" & snv$gene_symbol %in% pli$gene[which(pli$oe_mis_upper < 0.75)])), ]
# 
# ### read cnv data
# snv.osc <- snv[snv$Ethnicity == "EUR", ]#[snv$OSC1_OSC2 == "OSC1", ]
# snv.osc$freq_max <- 100*pmax(snv.osc$gnomAD_max, snv.osc$ALT_frq) #pmax(snv.osc$eur_ALT_frq, pmax(snv.osc$gnomAD_exome_AF, snv.osc$gnomAD_genome_AF, na.rm=T), na.rm=T)
# 
# snv.osc <- snv.osc[which(snv.osc$LOF == "LOF" | snv.osc$MSN %in% c("DMSN")), ]
# snv.osc <- snv.osc[which(is.na(snv.osc$SegDup)), ]
# snv.osc$varid <- paste(snv.osc$SID, snv.osc$chr, snv.osc$start, snv.osc$ref_allele, snv.osc$alt_allele, sep="#")
# snv.osc <- snv.osc[snv.osc$chr %in% paste0("chr", 1:22), ]
sampleinfo <- sampleinfo[!sampleinfo$TAG_id %in% c(outliers, samples_with_CNV3mb), ]

# cnv$SID[intersect(grep("KIF20B", cnv$Gene_Exons), which(cnv$max_freq < 0.5 & cnv$cnvFlag == "Gain")), ]
# cnv[intersect(grep("SAMD11", cnv$Gene_Exons), which(cnv$max_freq > 1 & cnv$max_freq < 5 & cnv$cnvFlag == "Gain")), ]
# 
# snv.osc[which(snv.osc$gene_symbol == "DST" & snv.osc$LOF == "LOF" & snv.osc$freq_max < 0.5), ]
# sum(snv.osc$gene_symbol == "GNAI1" & snv.osc$MSN == "DMSN" & snv.osc$freq_max < 0.5)
# snv.osc[which(snv.osc$gene_symbol == "CHD6" & snv.osc$MSN == "DMSN" & snv.osc$freq_max > 1 & snv.osc$freq_max < 5), ]
# kif20b <- cnv$SID[intersect(grep("KIF20B", cnv$Gene_Exons), which(cnv$max_freq < 0.5 & cnv$cnvFlag == "Gain")) ]
# boxplot(sampleinfo$itssrt_mf_tscores[sampleinfo$TAG_id %in% kif20b], sampleinfo$itssrt_mf_tscores[!sampleinfo$TAG_id %in% kif20b], names = c("KIF20B", "others"))

freq <- list(c(0, 0.5), c(1,5))
load("../../gsOSC.RData")
gsOSC <- gsOSC[-which(names(gsOSC) %in% c("ExprNov_BrainFeAd_sp", "Ilmn_BM.log2fpkm", "BspanML_lg2rpkm0.93", "BspanLA_lg2rpkm.MIN", "blueModule"))]

global.feats <- c("cnv_size_gain", "cnv_size_loss", "gene_count_gain", "gene_count_loss") #, "Total_LOF", "Total_DMSN")

for(freq_range in freq){
  ### filter CNVs
  cnv.tmp <- cnv[cnv$max_freq >= freq_range[1] & cnv$max_freq < freq_range[2], ]
  cnv.tmp <- cnv.tmp[cnv.tmp$SID %in% sampleinfo$TAG_id, ]
  cnvs <- CNVSet(cnv.tmp$SID, cnv.tmp$chr, cnv.tmp$start, cnv.tmp$end, tolower(cnv.tmp$cnvFlag))
  cnvs <- cnvs[cnvs$chr %in% paste0("chr", c(1:22)), ]
  cnvs$size <- cnvs$end - cnvs$start + 1
  
  ### filter SNVs
  # snv.tmp <- snv.osc[snv.osc$freq_max >= freq_range[1] & snv.osc$freq_max < freq_range[2], ]
  # snv.tmp$type <- ifelse(snv.tmp$LOF == "LOF" | snv.tmp$Splicing == "Splicing", "LOF", "DMSN")
  # snv.tmp$sample <- snv.tmp$SID
  
  dt.matrix <- getCNVGSMatrix(cnvs, genes, gsOSC) #, getSNVGSMatrix(snv.tmp, gsOSC), by = "sample", all = T)
  dt.matrix <- merge(dt.matrix, sampleinfo, by.x = "sample", by.y = "TAG_id", all.y = T)
  dt.matrix$gender_genetic <- factor(dt.matrix$gender_genetic, levels = c("Female", "Male"))
  
  for(i in 2:65){
    dt.matrix[is.na(dt.matrix[, i]), i] <- 0
  }
  
  
  covariates <- c("Study", "PCA1", "PCA2", "PCA3", "Plate")
  
  th.matrix <- dt.matrix
  
  global.out <- data.frame()
  for(global.feat in global.feats){
    ref.global <- glm(sprintf("gender_genetic ~ %s", paste(c(covariates), collapse = " + ")), th.matrix, family = binomial)
    add.global <- glm(sprintf("gender_genetic ~ %s", paste(c(covariates, global.feats), collapse = " + ")), th.matrix, family = binomial)
    conf <- confint(add.global)
    sm <- summary(add.global)
    
    global.out <- rbind(global.out, data.frame("freq" = paste(freq_range, collapse = "-"),
                                               "p" = sm$coefficients[global.feat, "Pr(>|z|)"], "coefficient" =  add.global$coefficients[global.feat],
                                               "coefficient.low" =  conf[global.feat, 1], "coefficient.up" =  conf[global.feat, 2],
                                               "feat" = global.feat))
  }
  
  dir.out <- sprintf("../results/EAS/%s/%s", paste(freq_range, collapse = "-"), Sys.Date())
  
  if(!dir.exists(dir.out)){
    dir.create(dir.out, recursive = T)
  }
  
  write.table(global.out, sprintf("%s/global_burden_cnv_sex_diff.tsv", dir.out), sep="\t", row.names=F, quote=F, col.names=T)
  
}
