library(pheatmap)

rmsys <- read.csv("RM_sys.csv", row.names=1)
rmsyssum <- read.csv("RM_sys_sum.csv", row.names=1)
defsys <- read.csv("defensesys.csv", row.names=1)
plasreps <- read.csv("summary_reps.tsv", row.names=1, header=TRUE, sep="\t", check.names=FALSE)
plasreps <- t(plasreps)
order <- c("PF_EC19", "PF_EC22", "PF_EC12", "PF_EC07", "PF_EC02", "PF_EC03", "PF_EC10", "PF_EC21", "PF_EC16", "PF_EC25",
                "PF_KPN10", "PF_KPN14", "PF_KPN01", "PF_KPN16", "PF_KQ01", "PF_KPN09", "PF_KPN12", "PF_KPN02", "PF_KPN11", "PF_KQ04")

rmsys <- rmsys[, order]
rmsyssum <- rmsyssum[, order]
defsys <- defsys[, order]
plasreps <- plasreps[, order]


pheatmap(rmsys, color=c("gray95", "orange"), cluster_rows=F, cluster_cols=F, border_color="white", angle_col=45, cellwidth=18, cellheight=14, legend=F)
pheatmap(rmsyssum, color=c("gray95", "orange", "orange3", "darkorange3", "darkorange4"), cluster_rows=F, cluster_cols=F, border_color="white", angle_col=45, cellwidth=18, cellheight=14, legend=F)
pheatmap(defsys, color=c("gray95", "mediumpurple1", "mediumpurple4"), cluster_rows=F, cluster_cols=F, border_color="white", angle_col=45, cellwidth=18, cellheight=14, legend=F)
pheatmap(plasreps, color=c("gray95", "darkolivegreen3"), cluster_rows=F, cluster_cols=F, border_color="white", angle_col=45, cellwidth=18, cellheight=14, legend=F)
