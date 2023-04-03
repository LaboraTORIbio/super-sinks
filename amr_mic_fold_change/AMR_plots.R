library(pheatmap)
library(RColorBrewer)

amr <- read.csv("summary_amr.tsv", row.names=1, header=TRUE, sep="\t", check.names=FALSE)
amr <- t(amr)

order_amc <- c("PF_EC25", "PF_EC10", "PF_EC23", "PF_EC24", "PF_EC07", "PF_EC11", "PF_EC20", "PF_EC21", "PF_EC22", "PF_EC06", "PF_EC02", "PF_EC03", "PF_EC16", "PF_EC05", "PF_EC14", "PF_KPN18", "PF_KPN16", "PF_KPN19", "PF_KPN09", "PF_KQ02", "PF_KQ03", "PF_KPN12", "PF_KPN01", "PF_KPN04", "PF_KPN03", "PF_KPN13", "PF_KPN11", "PF_KPN15", "PF_KPN20", "PF_KQ04", "PF_KPN06", "PF_KPN08", "PF_KV01")
order_ert <- c("PF_EC25", "PF_EC21", "PF_EC10", "PF_EC20", "PF_EC02", "PF_EC05", "PF_EC23", "PF_EC07", "PF_EC06", "PF_EC03", "PF_EC16", "PF_EC14", "PF_EC24", "PF_EC11", "PF_EC22", "PF_KPN16", "PF_KPN18", "PF_KPN04", "PF_KPN15", "PF_KPN09", "PF_KQ03", "PF_KPN12", "PF_KPN01", "PF_KPN03", "PF_KPN11", "PF_KPN20", "PF_KQ04", "PF_KPN19", "PF_KQ02", "PF_KV01", "PF_KPN13", "PF_KPN06", "PF_KPN08")
order_imp <- c("PF_EC20", "PF_EC25", "PF_EC10", "PF_EC05", "PF_EC07", "PF_EC06", "PF_EC03", "PF_EC14", "PF_EC24", "PF_EC21", "PF_EC02", "PF_EC23", "PF_EC16", "PF_EC11", "PF_EC22", "PF_KQ04", "PF_KQ03", "PF_KPN20", "PF_KPN19", "PF_KPN13", "PF_KPN12", "PF_KPN01", "PF_KQ02", "PF_KV01", "PF_KPN08", "PF_KPN09", "PF_KPN16", "PF_KPN18", "PF_KPN04", "PF_KPN11", "PF_KPN06", "PF_KPN15", "PF_KPN03")
order_mer <- c("PF_EC20", "PF_EC10", "PF_EC24", "PF_EC25", "PF_EC05", "PF_EC06", "PF_EC14", "PF_EC02", "PF_EC23", "PF_EC07", "PF_EC03", "PF_EC21", "PF_EC11", "PF_EC22", "PF_EC16", "PF_KPN08", "PF_KPN03", "PF_KQ03", "PF_KPN20", "PF_KPN19", "PF_KPN12", "PF_KPN01", "PF_KV01", "PF_KPN16", "PF_KPN18", "PF_KPN04", "PF_KPN15", "PF_KQ04", "PF_KPN11", "PF_KPN06", "PF_KPN13", "PF_KQ02", "PF_KPN09")

amr_amc <- amr[, order_amc]
amr_ert <- amr[, order_ert]
amr_imp <- amr[, order_imp]
amr_mer <- amr[, order_mer]

pheatmap(amr_amc,color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100), cluster_rows=F, cluster_cols=F, border_color="white", angle_col=45)
pheatmap(amr_ert,color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100), cluster_rows=F, cluster_cols=F, border_color="white", angle_col=45)
pheatmap(amr_imp,color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100), cluster_rows=F, cluster_cols=F, border_color="white", angle_col=45)
pheatmap(amr_mer,color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100), cluster_rows=F, cluster_cols=F, border_color="white", angle_col=45)