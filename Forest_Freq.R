### the first line only need to run once to install the GSBurden package
#install.packages(c("rstudioapi", "GenomicRanges", "BiasedUrn", "MASS", "ordinal"))

library(GenomicRanges)
library(GSBurden)
library(ggplot2)

#################################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sampleinfo <- readxl::read_excel("../phenotypes/Supplementary Table S1.v13.xlsx", sheet = 1, 
                                 col_types = c(rep("text", 11), rep("numeric", 29)))
sampleinfo$Plate <- paste0("Plate",gsub("Plate\\ ", "", sampleinfo$Plate))

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
cnv <- cnv[cnv$max_freq < 5, ]

samples_with_CNV3mb <- cnv$SID[cnv$size >= 3000000]
writeLines(samples_with_CNV3mb, "samples_with_CNV3mb.txt")
outliers <- readLines("../qcscripts/outliers.inAll.cnv.txt") #, readLines("../qcscripts/outliers.inAll.LoF.txt"))

remove <- c(outliers, samples_with_CNV3mb)

#### get numbers for supplementary mat
sum(sampleinfo$TAG_id[!is.na(sampleinfo$swan_mf_tscores)] %in% remove)
sum(sampleinfo$TAG_id[!is.na(sampleinfo$tocs_mf_tscores)] %in% remove)
sum(sampleinfo$TAG_id[!is.na(sampleinfo$itssrt_mf_tscores)] %in% remove)

sampleinfo <- sampleinfo[!sampleinfo$TAG_id %in% remove, ]

table(sampleinfo$ethnicity[!is.na(sampleinfo$swan_mf_tscores)])
table(sampleinfo$ethnicity[!is.na(sampleinfo$tocs_mf_tscores)])
table(sampleinfo$ethnicity[!is.na(sampleinfo$itssrt_mf_tscores)])

sum(!is.na(sampleinfo$swan_mf_tscores) & !is.na(sampleinfo$tocs_mf_tscores) & !is.na(sampleinfo$itssrt_mf_tscores))

sum(!is.na(sampleinfo$swan_mf_tscores) & !is.na(sampleinfo$tocs_mf_tscores) & !is.na(sampleinfo$itssrt_mf_tscores))/nrow(sampleinfo)
#####

sampleinfo <- sampleinfo[sampleinfo$ethnicity == "EUR", ]
### read snv data
# snv <- data.table::fread("../snvs/OSC1_OSC2.All.SNVs.hg19.Eur.EastAsian.10prct.rm_too_high_internal.txt", data.table = F)
# snv <- snv[which((snv$LOF == "LOF" & snv$gnomAD_oe_lof_upper < 0.35) |
#                    (snv$MSN == "DMSN" & snv$LOF != "LOF" & snv$gene_symbol %in% pli$gene[which(pli$oe_mis_upper < 0.75)])), ]
# 
# # snv <- snv[which(snv$gnomAD_pLI > 0.9), ]
# ### read cnv data
# snv.osc <- snv[snv$Ethnicity == "EUR", ]#[snv$OSC1_OSC2 == "OSC1", ]
# snv.osc$freq_max <- 100* snv.osc$gnomAD_max #pmax(snv.osc$eur_ALT_frq, pmax(snv.osc$gnomAD_exome_AF, snv.osc$gnomAD_genome_AF, na.rm=T), na.rm=T)
# 
# snv.osc <- snv.osc[which(snv.osc$LOF == "LOF" | snv.osc$MSN %in% c("DMSN")), ]
# snv.osc <- snv.osc[which(is.na(snv.osc$SegDup)), ]
# snv.osc$varid <- paste(snv.osc$SID, snv.osc$chr, snv.osc$start, snv.osc$ref_allele, snv.osc$alt_allele, sep="#")
# snv.osc <- snv.osc[snv.osc$chr %in% paste0("chr", 1:22), ]

freq <- list(c(0,0.1), c(0.1, 0.5), c(0.5,1), c(1,5))
load("../../gsOSC.RData")
gsOSC <- gsOSC[-which(names(gsOSC) %in% c("ExprNov_BrainFeAd_sp", "Ilmn_BM.log2fpkm", "BspanML_lg2rpkm0.93", "BspanLA_lg2rpkm.MIN", "blueModule"))]

dt.out <- data.frame()

feats <- c("gene_count_gain", "gene_count_loss") #, "Total_LOF", "Total_DMSN")
for(freq_range in freq){
  
  ### filter CNVs
  cnv.tmp <- cnv[cnv$max_freq >= freq_range[1] & cnv$max_freq < freq_range[2], ]
  cnv.tmp <- cnv.tmp[cnv.tmp$SID %in% sampleinfo$TAG_id, ]
  cnvs <- CNVSet(cnv.tmp$SID, cnv.tmp$chr, cnv.tmp$start, cnv.tmp$end, tolower(cnv.tmp$cnvFlag))
  cnvs <- cnvs[cnvs$chr %in% paste0("chr", c(1:22)), ]
  message(nrow(cnvs))
  ### filter SNVs
  # snv.tmp <- snv.osc[snv.osc$freq_max >= freq_range[1] & snv.osc$freq_max < freq_range[2], ]
  # snv.tmp$type <- ifelse(snv.tmp$LOF == "LOF" | snv.tmp$Splicing == "Splicing", "LOF", "DMSN")
  # snv.tmp$sample <- snv.tmp$SID
  
  dt.matrix <- getCNVGSMatrix(cnvs, genes, gsOSC) #, getSNVGSMatrix(snv.tmp, gsOSC), by = "sample", all = T)
  dt.matrix <- merge(dt.matrix, sampleinfo, by.x = "sample", by.y = "TAG_id", all.y = T)
  ### checking NA
  
  for(i in 2:65){
    dt.matrix[is.na(dt.matrix[, i]), i] <- 0
  }

  for(measure in c("swan_mf_tscores", "swan_ia_mf_tscores", "swan_hi_mf_tscores", 
                   "tocs_mf_tscores", "itssrt_mf_tscores", "crtsdt_mf_tscores")){
    
    covariates <- c("Study", "PCA1", "PCA2", "PCA3", "Plate")
    
    th.matrix <- dt.matrix[!is.na(dt.matrix[, measure]), ]
    
    if(measure %in% c("itssrt_mf_tscores", "crtsdt_mf_tscores")){
      th.matrix$stim_med <- th.matrix$"Stim..Meds"
      covariates <- c(covariates, "stim_med")
      th.matrix <- th.matrix[!is.na(th.matrix$"Stim..Meds"), ]
    }
    
    ref <- lm(sprintf("%s ~ %s", measure, paste(covariates, collapse = " + ")), th.matrix)
    add <- lm(sprintf("%s ~ %s + %s", measure, paste(covariates, collapse = " + "), paste(feats, collapse = " + ")), th.matrix)
    ano <- anova(add, ref, test = "Chisq")
    conf <- confint(add)
    
    global.test.out <- data.frame(measure, "freq" = paste(freq_range, collapse = "-"),
                                  "p" = ano$`Pr(>Chi)`[2], "coefficient" =  add$coefficients[feats],
                                  "coefficient.low" =  conf[feats, 1], "coefficient.up" =  conf[feats, 2], "global" = feats)

    dt.out <- rbind(dt.out, global.test.out)
  }
}

write.table(dt.out, "../results/EUR/freq_forest_plot_cnvs.tsv", sep="\t", row.names=F, quote=F, col.names=T)
dt.out$freq <- factor(dt.out$freq, levels = unique(dt.out$freq))

# global.map <- data.frame("global" = c("gene_count_loss", "gene_count_gain"), #, "Total_LOF", "Total_DMSN"),
#                          "plot_global" = c("Deletions", "Duplications")) #, "LoFs", "Damaging missenses"))
dt.out$plot_global <- ifelse(dt.out$global == "gene_count_loss", "Deletions", "Duplications")

ggplot(dt.out, aes(x = freq, y = coefficient, color = plot_global)) + geom_hline(yintercept = 0, lty = 2, lwd = .5) +
  geom_point(aes(fill = plot_global), color = "black", position = position_dodge(width = .5), shape = 21) +
  geom_errorbar(aes(ymin = coefficient.low, ymax = coefficient.up), position = position_dodge(width = .5), width = .2) + 
  coord_flip(ylim = c(-1, 1)) + facet_wrap(measure~.) + theme_bw() + theme(legend.position = "top") +
  scale_color_manual(values = c("Deletions"="red", "Duplications"="darkblue"))  + 
  scale_fill_manual(values = c("Deletions"="red", "Duplications"="darkblue")) + xlab("Frequency %")


ggsave("../results/EUR/freq_forest_plot_cnvs.png", width = 7, height = 5)

#####
dt.out$coeff.sign <- sign(dt.out$coefficient)
coeff <- data.frame()

for(freq in unique(dt.out$freq)){
  for(freq2 in unique(dt.out$freq)){
    match <- sum(dt.out$coeff.sign[dt.out$freq == freq] == dt.out$coeff.sign[dt.out$freq == freq2])
    denom <- match + (sum(dt.out$coeff.sign[dt.out$freq == freq] != dt.out$coeff.sign[dt.out$freq == freq2]))
    jacc <- signif(match/denom, digits = 2)
    jacc <- cor(dt.out$coefficient[dt.out$freq == freq], dt.out$coefficient[dt.out$freq == freq2], method = "spearman")
    coeff <- rbind(coeff, data.frame(freq, freq2, jacc))
  }
}
coeff
write.table(coeff, "../results/EUR/freq_forest_coeff.tsv", sep="\t", row.names=F, quote=F, col.names=T)
  