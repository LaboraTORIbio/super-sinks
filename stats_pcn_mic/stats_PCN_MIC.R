library("ggplot2")
library("ggpubr")

pf_tc_mic <- read.table("PF_TC_MICs.csv", header=TRUE)
tc_fc_pcn <- read.table("TC_MICs_FC_PCN.csv", header=TRUE)


### Comparing mean MICs (Fig. 1A)

# Comparing initial MIC (WT) E. coli vs Klebsiella

shapiro.test(pf_tc_mic[pf_tc_mic$Species == "E_coli" & pf_tc_mic$Type == "WT", ]$MIC_AMC)
shapiro.test(pf_tc_mic[pf_tc_mic$Species == "Klebsiella_spp" & pf_tc_mic$Type == "WT", ]$MIC_AMC)
wilcox.test(pf_tc_mic[pf_tc_mic$Species == "E_coli" & pf_tc_mic$Type == "WT", ]$MIC_AMC, pf_tc_mic[pf_tc_mic$Species == "Klebsiella_spp" & pf_tc_mic$Type == "WT", ]$MIC_AMC, paired = FALSE, alternative = "two.sided")

shapiro.test(pf_tc_mic[pf_tc_mic$Species == "E_coli" & pf_tc_mic$Type == "WT", ]$MIC_ERT)
shapiro.test(pf_tc_mic[pf_tc_mic$Species == "Klebsiella_spp" & pf_tc_mic$Type == "WT", ]$MIC_ERT)
wilcox.test(pf_tc_mic[pf_tc_mic$Species == "E_coli" & pf_tc_mic$Type == "WT", ]$MIC_ERT, pf_tc_mic[pf_tc_mic$Species == "Klebsiella_spp" & pf_tc_mic$Type == "WT", ]$MIC_ERT, paired = FALSE, alternative = "two.sided")

shapiro.test(pf_tc_mic[pf_tc_mic$Species == "E_coli" & pf_tc_mic$Type == "WT", ]$MIC_IMP)
shapiro.test(pf_tc_mic[pf_tc_mic$Species == "Klebsiella_spp" & pf_tc_mic$Type == "WT", ]$MIC_IMP)
wilcox.test(pf_tc_mic[pf_tc_mic$Species == "E_coli" & pf_tc_mic$Type == "WT", ]$MIC_IMP, pf_tc_mic[pf_tc_mic$Species == "Klebsiella_spp" & pf_tc_mic$Type == "WT", ]$MIC_IMP, paired = FALSE, alternative = "two.sided")

shapiro.test(pf_tc_mic[pf_tc_mic$Species == "E_coli" & pf_tc_mic$Type == "WT", ]$MIC_MER)
shapiro.test(pf_tc_mic[pf_tc_mic$Species == "Klebsiella_spp" & pf_tc_mic$Type == "WT", ]$MIC_MER)
wilcox.test(pf_tc_mic[pf_tc_mic$Species == "E_coli" & pf_tc_mic$Type == "WT", ]$MIC_MER, pf_tc_mic[pf_tc_mic$Species == "Klebsiella_spp" & pf_tc_mic$Type == "WT", ]$MIC_MER, paired = FALSE, alternative = "two.sided")

# Comparing final MIC (TC) E. coli vs Klebsiella

shapiro.test(tc_fc_pcn[tc_fc_pcn$Species == "E_coli", ]$MIC_AMC)
shapiro.test(tc_fc_pcn[tc_fc_pcn$Species == "Klebsiella_spp", ]$MIC_AMC)
wilcox.test(MIC_AMC ~ Species, data = tc_fc_pcn, paired = FALSE, alternative = "two.sided")

shapiro.test(tc_fc_pcn[tc_fc_pcn$Species == "E_coli", ]$MIC_ERT)
shapiro.test(tc_fc_pcn[tc_fc_pcn$Species == "Klebsiella_spp", ]$MIC_ERT)
wilcox.test(MIC_ERT ~ Species, data = tc_fc_pcn, paired = FALSE, alternative = "two.sided")

shapiro.test(tc_fc_pcn[tc_fc_pcn$Species == "E_coli", ]$MIC_IMP)
shapiro.test(tc_fc_pcn[tc_fc_pcn$Species == "Klebsiella_spp", ]$MIC_IMP)
wilcox.test(MIC_IMP ~ Species, data = tc_fc_pcn, paired = FALSE, alternative = "two.sided")

shapiro.test(tc_fc_pcn[tc_fc_pcn$Species == "E_coli", ]$MIC_MER)
shapiro.test(tc_fc_pcn[tc_fc_pcn$Species == "Klebsiella_spp", ]$MIC_MER)
wilcox.test(MIC_MER ~ Species, data = tc_fc_pcn, paired = FALSE, alternative = "two.sided")


### Correlation of PCN and final MIC of AMC

cor.test(tc_fc_pcn$MIC_AMC, tc_fc_pcn$plas.chr123_median, method = "spearman")

ggscatter(tc_fc_pcn, x = "MIC_AMC", y = "plas.chr123_median", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "MIC AMC", ylab = "pOXA-48 copy number",
          add.params = list(color = "black", fill = "lightgray"),
          size = 4, cor.coef.size = 6, color="Species", alpha=0.5)
