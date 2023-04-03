
library(tidyverse)
library(ggpubr)
library(readxl)

#Figure 1A. MICs

MICS<- read_excel("source_data.xlsx", 
                        sheet = "MICs (Fig. 1A)")
MICS %>%
  ggplot(aes(fill=Species, reorder(Type, MIC_AMC), MIC_AMC))+ # Change here the antibiotic to plot
  geom_line(aes(group=Strain), color="grey50")+
  scale_fill_manual(values=c("#f4646dff", "#009fe3ff"))+
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1, dotsize=.88, binwidth =.7)+
  scale_y_continuous(trans= ('log2'), breaks =c(16,128,1024,8192,65536))+
  facet_grid(~Species)+
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank())

#Figure 2A. One recipient conjugation frequencies

CC_single <- read_excel("source_data.xlsx", 
                        sheet = "Single conjugations (Fig. 2A)")

ggplot(CC_single, aes(x=Recipient, y=Value, color=Species, shape=Replicate)) +
  geom_jitter(position=position_dodge(0.8))+
  scale_color_manual(values=c("#f4646dff", "#009fe3ff"))+
  scale_shape_identity()+
  scale_y_continuous(trans = 'log10', limits = c(1e-08, 1e-00), 
                     breaks = c(1e-08, 1e-06, 1e-04, 1e-02, 1e-00))+
  facet_grid(~Donor)+
  theme_bw()

#Figure 2B. Recipient pairs conjugation frequencies

CC_pairs <- read_excel("source_data.xlsx", 
                       sheet = "Paired conjugations (Fig. 2B)")

ggplot(CC_pairs, aes(x=Pair_recipient, y=Value, color=Species, shape=Replicate)) +
  geom_jitter(position=position_dodge(0.8))+
  scale_color_manual(values=c("#f4646dff", "#009fe3ff"))+
  scale_shape_identity()+
  scale_y_continuous(trans = 'log10', limits = c(1e-08, 1e-00), 
                     breaks = c(1e-08, 1e-06, 1e-04, 1e-02, 1e-00))+
  facet_grid(~Donor)+
  theme_bw()

#Figure 2C. Recipients pool conjugation frequencies

CC_pool <- read_excel("source_data.xlsx", 
                      sheet = "Pool conjugations (Fig. 2C)")

ggplot(CC_pool, aes(x=Replicate, y=Value, color=Species_recipient, shape=Replicate)) +
  geom_jitter(position=position_dodge(0.8))+
  scale_color_manual(values=c("#f4646dff", "#009fe3ff"))+
  scale_shape_identity()+
  scale_y_continuous(trans = 'log10', limits = c(1e-08, 1e-00), 
                     breaks = c(1e-08, 1e-06, 1e-04, 1e-02, 1e-00))+
  scale_x_discrete('Replicate')+
  facet_grid(~Donor)+
  theme_bw()

#Figure 3A. Mouse gut colonization levels

cfu <- read_excel("source_data.xlsx", 
                         sheet = "Colonization levels (Fig. 3A)")

ggplot(cfu, aes(x=Recipient, y=logCFU, color=Species, shape=Replicate)) +
  geom_jitter(position=position_dodge(0.8))+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5, color= "black")+
  scale_color_manual(values=c("#f4646dff", "#009fe3ff"))+
  scale_shape_identity()+
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2, 4, 6, 8,10))+
  facet_grid(Sample ~ Type)+
  theme_bw()

#Figure 3B. In vivo conjugation frequencies

conj_mice <- read_excel("source_data.xlsx", 
                        sheet = "Conjugation freqs (Fig. 3B)")

ggplot(conj_mice, aes(x=Recipient, y=Value, color=Species, shape=Replicate)) +
  geom_jitter(position=position_dodge(0.8))+
  #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #geom = "crossbar", width = 0.5, color= "black")+
  scale_color_manual(values=c("#f4646dff", "#009fe3ff"))+
  scale_shape_identity()+
  scale_y_continuous(trans = 'log10', limits = c(1e-08, 1e-00), 
                     breaks = c(1e-08, 1e-06, 1e-04, 1e-02, 1e-00))+
  facet_grid(~Sample)+
  theme_bw()


#Statistical analysis

#Paired Wilcoxon test for recipient pairs conjugations

pairwise.wilcox.test(CC_pairs$Value, CC_pairs$Species , paired = TRUE, alternative= "greater", p.adjust.method = "holm")

#Kruskal-wallis test for single and pool conjugations in vitro and in vivo conjugations

kruskal.test(Value ~ interaction(Donor,Species_recipient), data = CC_pool %>% filter(Donor=="beta3914"))

