library(pheatmap)
library(RColorBrewer)
library(tidyr)

order <- c("EC22", "EC12", "EC07", "EC02", "EC19", "EC03", "EC10", "EC21", "EC16", "EC25",
           "KPN10", "KPN01", "KPN14", "KPN16", "KQ01", "KPN09", "KPN12", "KQ04", "KPN02", "KPN11")

caps_sets <- read.table("caps_sets/wGRR_genomes_Twoway.txt")
mat <- pivot_wider(caps_sets, names_from = V2, values_from = V3)
mat <- as.matrix(mat)
matrix <- mat[,-1]
rownames(matrix) <- mat[,1]
class(matrix) <- "numeric"
mator <- matrix[order, order]
colnames(mator) <- c("PF_EC22", "PF_EC12", "PF_EC07", "PF_EC02", "PF_EC19", "PF_EC03", "PF_EC10", "PF_EC21", "PF_EC16", "PF_EC25",
                     "PF_KPN10", "PF_KPN01", "PF_KPN14", "PF_KPN16", "PF_KQ01", "PF_KPN09", "PF_KPN12", "PF_KQ04", "PF_KPN02", "PF_KPN11")
rownames(mator) <- c("PF_EC22", "PF_EC12", "PF_EC07", "PF_EC02", "PF_EC19", "PF_EC03", "PF_EC10", "PF_EC21", "PF_EC16", "PF_EC25",
                     "PF_KPN10", "PF_KPN01", "PF_KPN14", "PF_KPN16", "PF_KQ01", "PF_KPN09", "PF_KPN12", "PF_KQ04", "PF_KPN02", "PF_KPN11")



## Histograms and heatmaps
colors <- colorRampPalette(brewer.pal(9, "Blues"))(140)
pheatmap(mator, cluster_rows=F, cluster_cols=F, color=colors, border_color="white", angle_col=0, breaks=seq(0,1, length.out=100), cellwidth=16, cellheight=14.8)
