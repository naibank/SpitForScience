library(data.table)
library(ggplot2)

### read cnv data of the first batch of the cohort
cnv.osc1 <- read.delim("../cnvs/Fragmentation_Merge.cnv.annotation20190807.Annot.txt", stringsAsFactors = F)

### retain stringent call sets
cnv.osc1 <- cnv.osc1[cnv.osc1$stringentFlag == "stringent", ]
cnv.osc1 <- subset(cnv.osc1, size >= 10000 & probeNum >= 5 &
                     unclean < 70 & segdup < 70 &  repeats < 70 &
                     DGV.OnePrcnt.Clusters == "No.Cluster(s)=0") 
cnv.osc1$Against.East.Asian.Int.Freq <- cnv.osc1$Against.East.Asian.Int.Freq - min(cnv.osc1$Against.East.Asian.Int.Freq, na.rm = T)
cnv.osc1$Against.Eur.Int.Freq <- cnv.osc1$Against.Eur.Int.Freq - min(cnv.osc1$Against.Eur.Int.Freq, na.rm = T)

### calculate max frequency from EUR and East Asian subsets, internally and externally
cnv.osc1$max_freq <- pmax(cnv.osc1$Against.All.TCAG.Conts.Freq, pmax(cnv.osc1$Against.Eur.TCAG.Cont.Freq, 
                                                                     pmax(cnv.osc1$Against.Chinese.Conts.Freq,
                                                                          pmax(cnv.osc1$Against.Eur.Int.Freq, cnv.osc1$Against.East.Asian.Int.Freq, na.rm = T), na.rm = T), na.rm = T), na.rm = T)

### read cnv data of the second batch of the cohort
cnv.osc2 <- read.delim("../cnvs/sfs3141passQC.TAGID_gender.cnv_annotation20190807.v2.txt", stringsAsFactors = F)
cnv.osc2 <- cnv.osc2[cnv.osc2$stringentFlag == "stringent", ]
cnv.osc2 <- subset(cnv.osc2, size >= 10000 & probeNum >= 5 &
                     unclean < 70 & segdup < 70 &  repeats < 70 &
                     DGV.OnePrcnt.Clusters == "No.Cluster(s)=0")
cnv.osc2$Internal.Freq <- cnv.osc2$Internal.Freq-min(cnv.osc2$Internal.Freq, na.rm = T)

### calculate max frequency from EUR and East Asian subsets, internally and externally
cnv.osc2$max_freq <- pmax(cnv.osc2$Against.TCAG.Conts.All.Freq, cnv.osc2$Internal.Freq, na.rm = T)

feat <- intersect(names(cnv.osc1), names(cnv.osc2))
cnv <- rbind(cnv.osc1[, feat], cnv.osc2[, feat])
cnv$max_freq[is.na(cnv$max_freq)] <- 0
cnv <- cnv[cnv$max_freq < 5, ]

# cnv <- cnv[cnv$Gene_Exons != "", ]
# all <- data.frame(table(snv$SID[!duplicated(snv$varid)]))
del <- cnv[cnv$cnvFlag == "Loss", ]
del <- data.frame(table(del$SID))
dup <- cnv[cnv$cnvFlag == "Gain", ]
dup <- data.frame(table(dup$SID))

# dt <- merge(all, lof, by = "Var1", all = T)
dt <- merge(del, dup, by = "Var1", all = T)
dt[is.na(dt)] <- 0
names(dt) <- c("SID", "DEL", "DUP")

write.table(dt, "cnv_var_count.tsv", sep="\t", row.names=F, quote=F, col.names=T)

dt <- read.delim("cnv_var_count.tsv", stringsAsFactors = F)
saminfo <- readxl::read_excel("../phenotypes/Phenotypes for S1 and S2 CNV Analysis-Jan 19 2022.xlsx", sheet = 2)
saminfo <- merge(saminfo, dt, by.x = "id", by.y = "SID", all.x = T)
# saminfo$All[is.na(saminfo$All)] <- 0
saminfo$DEL[is.na(saminfo$DEL)] <- 0
saminfo$DUP[is.na(saminfo$DUP)] <- 0
outliers <- c()
for(stdu in unique(saminfo$Study)){
  for(pop in unique(saminfo$ethn_white)){
    for(var in c("DEL", "DUP")){
      count <- saminfo[saminfo$Study == stdu & saminfo$ethn_white == pop, ]
      count$anscombe <- 2*sqrt(count[, var] + 3/8)
      cutoff <- mean(count$anscombe) + 3 * sd(count$anscombe)
      outliers <- union(outliers, count$id[count$anscombe > cutoff])
      message(sprintf("%s %s %s %s", stdu, pop, var, length(outliers)))
    }
  }
}
saminfo$outlier <- saminfo$id %in% outliers

# all <- ggplot(saminfo, aes(x = paste(sprintf("OSC%s",Study), ifelse(ethn_white == 1, "EUR", "EAS")), y = All)) + 
#   geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(fill = outlier), alpha = 2, shape = 21, color = "black", show.legend = F) +
#   xlab("") + theme_bw() + scale_fill_manual(values = c("grey", "red"))

del <- ggplot(saminfo, aes(x = paste(sprintf("Spit %s",Study), ifelse(ethn_white == 1, "EUR", "EAS")), y = DEL)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(aes(fill = outlier), width =.2, alpha = 2, size = 2, shape = 21, color = "black", show.legend = F) +
  xlab("") + theme_bw() + scale_fill_manual(values = c("grey", "red")) + ylab("Number of deletions (<5%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dup <- ggplot(saminfo, aes(x = paste(sprintf("Spit %s",Study), ifelse(ethn_white == 1, "EUR", "EAS")), y = DUP)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(aes(fill = outlier), width =.2, alpha = 2, size = 2, shape = 21, color = "black", show.legend = F) +
  xlab("") + theme_bw() + scale_fill_manual(values = c("grey", "red")) + ylab("Number of duplications (<5%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cowplot::plot_grid(del, dup, nrow = 1)
ggsave("CNV.QC.png", width = 6, height = 5)

writeLines(outliers, "outliers.inAll.cnv.txt")

sum(saminfo$id[!is.na(saminfo$swan_mf_tscores)] %in% outliers) #51
sum(saminfo$id[!is.na(saminfo$itssrt_mf_tscores)] %in% outliers) #46
sum(saminfo$id[!is.na(saminfo$tocs_mf_tscores)] %in% outliers) #51

