### the first line only need to run once to install the GSBurden package
#install.packages(c("rstudioapi", "GenomicRanges", "BiasedUrn", "MASS", "ordinal"))

library(GenomicRanges)
library(GSBurden)
library(ggplot2)

#################################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sampleinfo <- readxl::read_excel("../phenotypes/Supplementary Table S1.v13.xlsx", sheet = 1, 
                                 col_types = c(rep("text", 11), rep("numeric", 29)))
sampleinfo <- sampleinfo[sampleinfo$ethnicity == "EUR", ]
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

write.table(cnv, "stringent_CNVs_Spit12.tsv", sep="\t", row.names=F, quote=F, col.names=T)

samples_with_CNV3mb <- cnv$SID[cnv$size >= 3000000]
# writeLines(samples_with_CNV3mb, "samples_with_CNV3mb.txt")

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
  
  dir.out <- sprintf("../results/EUR/%s/%s", paste(freq_range, collapse = "-"), Sys.Date())
  
  if(!dir.exists(dir.out)){
    dir.create(dir.out, recursive = T)
  }
  
  write.table(global.out, sprintf("%s/global_burden_cnv_sex_diff.tsv", dir.out), sep="\t", row.names=F, quote=F, col.names=T)
  
}

eur.out <- data.frame()
for(freq_range in list(c(0, 0.5), c(1,5))){
  dir.in <- sprintf("../results/EUR/%s/%s", paste(freq_range, collapse = "-"), "2023-04-24")
  dt.tmp <- read.delim(sprintf("%s/global_burden_cnv_sex_diff.tsv", dir.in), stringsAsFactors = F)
  dt.tmp <- dt.tmp[dt.tmp$feat %in% global.feats, ]
  analysis <- paste("EUR", ifelse(freq_range == c(0, 0.5), "Rare", "Common"))
  dt.tmp$analysis <- analysis
  eur.out <- rbind(eur.out, dt.tmp)
}

eas.out <- data.frame()
for(freq_range in list(c(0, 0.5), c(1,5))){
  dir.in <- sprintf("../results/EUR/%s/%s", paste(freq_range, collapse = "-"), "2023-04-24")
  dt.tmp <- read.delim(sprintf("%s/global_burden_cnv_sex_diff.tsv", dir.in), stringsAsFactors = F)
  dt.tmp <- dt.tmp[dt.tmp$feat %in% global.feats, ]
  analysis <- paste("EAS", ifelse(freq_range == c(0, 0.5), "Rare", "Common"))
  dt.tmp$analysis <- analysis
  eas.out <- rbind(eas.out, dt.tmp)
}

dt.out <- rbind(eur.out, eas.out)
dt.out$coefficient <- signif(dt.out$coefficient, digits = 2)
dt.out$p <- signif(dt.out$p, digits = 2)
dt.out$coefficient.low <- signif(dt.out$coefficient.low, digits = 2)
dt.out$coefficient.up <- signif(dt.out$coefficient.up, digits = 2)
names(dt.out)[6] <- "variable"
write.table(dt.out, "../dataToShare/supplementary_sex_related_analysis.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# #### generate plots
# for(freq_range in freq){
#   
#   dir.in <- sprintf("../results/EUR/%s/%s", paste(freq_range, collapse = "-"), "2022-05-17")
#   dt.out <- read.delim(sprintf("%s/geneset_burden_cnv.tsv", dir.in), stringsAsFactors = F)
#   
#   qqnorm(dt.out$p); qqline(dt.out$p)
# 
#   for(measure in unique(dt.out$measure)){
#     ggplot(dt.out[dt.out$measure == measure, ], aes(x = geneset, y = coefficient, fill = BHFDR < 0.20 & p < 0.05)) +
#       geom_bar(stat = "identity", color = "black", width = .5) + theme_bw() + #coord_cartesian(ylim=c(-1, 5)) +
#       geom_errorbar(aes(ymin = coefficient.low, ymax = coefficient.up), size = .5, width = .2) +
#       geom_hline(yintercept = 0) + 
#       theme(panel.grid = element_blank(),
#             axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
#       scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "lightblue")) +
#       facet_wrap(type ~., nrow = 4, strip.position = "right", scales = "free_y") + ggtitle(measure)
#     
#     ggsave(sprintf("%s/%s_geneset_cnv.png", dir.in, measure), width = 12, height = 7)
#     
#   }
#   
#   ## merge plot
#   dt.out$sig <- ifelse(dt.out$BHFDR < 0.20 & dt.out$p < 0.05, "*", "")
#   
#   type_map <- data.frame("type" = c("DMSN", "gain", "LOF", "loss"), "plot_type" = c("Damaging missenses", "Duplications", "LoFs", "Deletions"))
#   score_map <- data.frame("measure" = c("crtsdt_mf_tscores", "itssrt_mf_tscores", "swan_hi_mf_tscores", "swan_ia_mf_tscores", "swan_mf_tscores", "tocs_mf_tscores"), 
#                           "plot_measure" = c("Response time variable", "SSRT", "Hyperactivity", "Inattention", "ADHD", "OCD"))
#   dt.out <- merge(dt.out, type_map, by = "type", all.x = T)
#   dt.out <- merge(dt.out, score_map, by = "measure")
#   dt.out$plot_type <- factor(dt.out$plot_type, levels = c("Deletions", "Duplications")) #, "LoFs", "Damaging missenses"))
#   dt.out$plot_measure <- factor(dt.out$plot_measure, levels = c("ADHD", "SSRT", "OCD", "Hyperactivity", "Inattention", "Response time variable"))
#   
#   sigset <- unique(dt.out$geneset[which(dt.out$sig == "*")])
#   main <- dt.out[dt.out$geneset %in% sigset & dt.out$plot_measure %in% c("SSRT", "ADHD", "OCD"), ]
#   
#   m <- ggplot(main, aes(x = geneset, y = coefficient, fill = plot_measure)) + geom_bar(stat = "identity", position = position_dodge(width = .75), color = "black", width = .75) + 
#     geom_errorbar(aes(ymin = coefficient.low, ymax = coefficient.up), position = position_dodge(width = .75), size = .5, width = .2) +
#     facet_wrap(plot_type~., nrow = 6, scales="free_y") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#                                                                               legend.title = element_blank(), legend.position = "top") + 
#     geom_text(aes(y = coefficient.up + 0.1, label = sig), size = 5, position = position_dodge(width = .75)) +geom_hline(yintercept = 0, lty = 2) +
#     scale_fill_manual(values = c("ADHD" = "red", "SSRT" = "lightblue", "OCD" = "yellow")) + ggtitle("Main phenotypes")
#  
#   sub <- dt.out[dt.out$geneset %in% sigset & dt.out$plot_measure %in% c("Hyperactivity", "Inattention", "Response time variable"), ]
#   
#   s <- ggplot(sub, aes(x = geneset, y = coefficient, fill = plot_measure)) + geom_bar(stat = "identity", position = position_dodge(width = .75), color = "black", width = .75) + 
#     geom_errorbar(aes(ymin = coefficient.low, ymax = coefficient.up), position = position_dodge(width = .75), size = .5, width = .2) +
#     facet_wrap(plot_type~., nrow = 6, scales="free_y") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#                                                                               legend.title = element_blank(), legend.position = "top") + 
#     geom_text(aes(y = coefficient.up + 0.1, label = sig), size = 5, position = position_dodge(width = .75)) +geom_hline(yintercept = 0, lty = 2) +
#     scale_fill_manual(values = c("Hyperactivity" = "salmon", "Inattention" = "indianred", "Response time variable" = "blue")) + ggtitle("Sub phenotypes")
#   
#   cowplot::plot_grid(m, s, nrow = 1)
#   ggsave(sprintf("%s/combined_geneset_cnv.png", dir.in), width = ifelse(length(sigset)*2 < 10, 10, length(sigset)*2), height = 9)
# }
# 
# cytoband <- read.delim("../../../PGC/cross_disorder_burden/cytoBand.txt", stringsAsFactors = F, header = F)
# cytoband.g <- GRanges(cytoband$V1, IRanges(cytoband$V2, cytoband$V3), "*")
# gene.g <- GRanges(genes$chr, IRanges(genes$start, genes$end), "*")
# olap <- data.frame(findOverlaps(gene.g, cytoband.g))
# olap$cytoband <- cytoband$V4[olap$subjectHits]
# olap <- aggregate(cytoband ~ queryHits, olap, paste, collapse = ";")
# genes$cytoband <- ""
# genes$cytoband[olap$queryHits] <- olap$cytoband
# 
# 
# #### loci test
# for(freq_range in freq){
#   ### filter CNVs
#   cnv.tmp <- cnv[cnv$max_freq >= freq_range[1] & cnv$max_freq < freq_range[2], ]
#   cnv.tmp <- cnv.tmp[cnv.tmp$SID %in% sampleinfo$TAG_id, ]
#   cnvs <- CNVSet(cnv.tmp$SID, cnv.tmp$chr, cnv.tmp$start, cnv.tmp$end, tolower(cnv.tmp$cnvFlag))
#   cnvs <- cnvs[cnvs$chr %in% paste0("chr", c(1:22)), ]
#   cnvs$size <- cnvs$end - cnvs$start + 1
#   cnvs.g <- GRanges(cnvs$chr, IRanges(cnvs$start, cnvs$end), "*")
#   
#   ### filter SNVs
#   # snv.tmp <- snv.osc[snv.osc$freq_max >= freq_range[1] & snv.osc$freq_max < freq_range[2], ]
#   # snv.tmp$type <- ifelse(snv.tmp$LOF == "LOF" | snv.tmp$Splicing == "Splicing", "LOF", "DMSN")
#   # snv.tmp$sample <- snv.tmp$SID
#   
#   dir.in <- sprintf("../results/EUR/%s/%s", paste(freq_range, collapse = "-"), "2022-05-20")
#   gs <- read.delim(sprintf("%s/geneset_burden_cnv.tsv", dir.in), stringsAsFactors = F)
#   gs <- gs[which(gs$BHFDR < 0.20 & gs$p < 0.05), ]
#   
#   measures <- unique(gs$measure)
#   for(measure in measures){
#     covariates <- c("Study", "PCA1", "PCA2", "PCA3", "Plate")
#     
#     th.samples <- data.frame(sampleinfo[!is.na(sampleinfo[, measure]), ])
# 
#     if(measure %in% c("itssrt_mf_tscores", "crtsdt_mf_tscores")){
#       th.samples$stim_med <- th.samples$"Stim..Meds"
#       covariates <- c(covariates, "stim_med")
#       th.samples <- th.samples[!is.na(th.samples$"Stim..Meds"), ]
#     }
#     
#     genesets <- unique(gs$geneset[gs$measure == measure])
#     th.genes <- genes[genes$enzid %in% unlist(gsOSC[genesets]), ]
#     th.genes.g <- GRanges(th.genes$chr, IRanges(th.genes$start, th.genes$end), "*")
#     olap <- data.frame(findOverlaps(cnvs.g, th.genes.g))
#     olap$entrez_id <- th.genes$enzid[olap$subjectHits]
#     olap$gene_symbol <- th.genes$gsymbol[olap$subjectHits]
#     olap$sample <- cnvs$sample[olap$queryHits]
#     olap$type <- cnvs$type[olap$queryHits]
#     olap$chr <- cnvs$chr[olap$queryHits]
#     olap$cytoband <- th.genes$cytoband[olap$subjectHits]
# 
#     # th.snvs <- snv.tmp[snv.tmp$entrez_id %in% unlist(gsOSC[genesets]), ]
#     # test.genes <- unique(rbind(olap[, c("gene_symbol", "entrez_id")], th.snvs[, c("gene_symbol", "entrez_id")]))
#     test.genes <- unique(olap[, c("gene_symbol", "entrez_id", "chr", "cytoband")])
#     test.genes[, c("genesets", "gain", "loss", # "LOF", "DMSN", 
#                    "gain.low", "loss.low", # "LOF.low", "DMSN.low", 
#                    "gain.up", "loss.up", #"LOF.up", "DMSN.up", 
#                    "pvalue", "loss.samples", "gain.samples")] <- NA
#     
#     for(i in 1:nrow(test.genes)){
#       test.genes$genesets[i] <- paste(genesets[which(sapply(sapply(gsOSC[genesets], "%in%", test.genes$entrez_id[i]), sum) > 0)], collapse=";")
#       
#       # th.samples$LOF <- as.numeric(th.samples$TAG_id %in% th.snvs$SID[th.snvs$entrez_id == test.genes$entrez_id[i] & th.snvs$LOF == "LOF"])
#       # th.samples$DMSN <- as.numeric(th.samples$TAG_id %in% th.snvs$SID[th.snvs$entrez_id == test.genes$entrez_id[i] & th.snvs$LOF != "LOF" &  th.snvs$MSN == "DMSN"])
#       th.samples$loss <- as.numeric(th.samples$TAG_id %in% olap$sample[olap$entrez_id == test.genes$entrez_id[i] & olap$type == "loss"])
#       th.samples$gain <- as.numeric(th.samples$TAG_id %in% olap$sample[olap$entrez_id == test.genes$entrez_id[i] & olap$type == "gain"])
#       
#       ref <- lm(sprintf("%s ~ %s", measure, paste(c(covariates), collapse = " + ")), th.samples)
#       add <- lm(sprintf("%s ~ %s", measure, paste(c(covariates, "gain", "loss"), collapse = " + ")), th.samples) #, "LOF", "DMSN"), 
# 
#       ano <- anova(add, ref, test = "Chisq")
#       conf <- confint(add)
#       
#       test.genes[i, c("gain", "loss")] <- add$coefficients[c("gain", "loss")]
#       test.genes[i, c("gain.low", "loss.low")] <- conf[c("gain", "loss"), 1]
#       test.genes[i, c("gain.up", "loss.up")] <- conf[c("gain", "loss"), 2]
#       
#       loss.samples <- paste(sort(th.samples$TAG_id[th.samples$loss == 1]), collapse = ";")
#       gain.samples <- paste(sort(th.samples$TAG_id[th.samples$gain == 1]), collapse = ";")
#       
#       test.genes[i, c("loss.samples")] <- loss.samples
#       test.genes[i, c("gain.samples")] <- gain.samples
#       
#       test.genes$pvalue[i] <- ano$`Pr(>Chi)`[2]
#     }
#     
#     test.genes <- test.genes[order(test.genes$pvalue), ]
#     test.genes$BHFDR <- p.adjust(test.genes$pvalue, method = "BH")
#     write.table(test.genes, sprintf("%s/%s_loci_test_cnv.tsv", dir.in, measure), sep="\t", row.names=F, quote=F, col.names=T)
#   }
# }
