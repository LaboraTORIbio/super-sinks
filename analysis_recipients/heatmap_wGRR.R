library(pheatmap)
library(RColorBrewer)
library(tidyr)

order <- c("EC02", "EC22", "EC07", "EC19", "EC12", "EC03", "EC10", "EC21", "EC16", "EC25",
           "KPN10", "KPN01", "KPN16", "KQ01", "KPN14", "KPN12", "KPN09", "KPN11", "KPN02", "KQ04")

caps_sets <- read.table("wGRR_genomes_Twoway.txt")
mat <- pivot_wider(caps_sets, names_from = V2, values_from = V3)
mat <- as.matrix(mat)
matrix <- mat[,-1]
rownames(matrix) <- mat[,1]
class(matrix) <- "numeric"
mator <- matrix[order, order]
colnames(mator) <- c("PF_EC02", "PF_EC22", "PF_EC07", "PF_EC19", "PF_EC12", "PF_EC03", "PF_EC10", "PF_EC21", "PF_EC16", "PF_EC25",
                     "PF_KPN10", "PF_KPN01", "PF_KPN16", "PF_KQ01", "PF_KPN14", "PF_KPN12", "PF_KPN09", "PF_KPN11", "PF_KPN02", "PF_KQ04")
rownames(mator) <- c("PF_EC02", "PF_EC22", "PF_EC07", "PF_EC19", "PF_EC12", "PF_EC03", "PF_EC10", "PF_EC21", "PF_EC16", "PF_EC25",
                     "PF_KPN10", "PF_KPN01", "PF_KPN16", "PF_KQ01", "PF_KPN14", "PF_KPN12", "PF_KPN09", "PF_KPN11", "PF_KPN02", "PF_KQ04")



## Histograms and heatmaps
colors <- colorRampPalette(brewer.pal(9, "Blues"))(140)
pheatmap(mator, cluster_rows=F, cluster_cols=F, color=colors, border_color="white", angle_col=0, breaks=seq(0,1, length.out=100), cellwidth=16, cellheight=14.8)
