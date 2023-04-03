library(ape)
library(phylolm)
library(phytools)


# Reading the tree files
tree_refseq <- read.tree(file = "../hmm_gene_markers/refseq_alignments/concat.mafft.trim.msa.treefile")
tree_refseq_ec <- read.tree(file = "../hmm_gene_markers/refseq_alignments_Ec/concat.mafft.trim.msa.treefile")

# Rooting trees
tree_refseq <- root(tree_refseq, "GCF_003204095.1", node = NULL, resolve.root = TRUE)
tree_refseq_ec <- root(tree_refseq_ec, "GCF_003204095.1", node = NULL, resolve.root = TRUE)

# Removing root branch
tree_refseq <- drop.tip(tree_refseq, "GCF_003204095.1")
tree_refseq_ec <- drop.tip(tree_refseq_ec, "GCF_003204095.1")

# Filling the null branches with a value >0
tree_refseq$edge.length <- replace(tree_refseq$edge.length, tree_refseq$edge.length==0, 0.000000000001)
tree_refseq_ec$edge.length <- replace(tree_refseq_ec$edge.length, tree_refseq_ec$edge.length==0, 0.000000000001)

# Rescaling tree
tree_refseq$edge.length <- tree_refseq$edge.length*(1/max(nodeHeights(tree_refseq)))
tree_refseq_ec$edge.length <- tree_refseq_ec$edge.length*(1/max(nodeHeights(tree_refseq_ec)))

# Exporting the tree nodes in order
write(tree_refseq$tip.label, "./treenodes_all.txt")
write(tree_refseq_ec$tip.label, "./treenodes_Ec.txt")


# From the terminal, run the script get_traits_table.py to generate the traits data frame
# Remember to have the lists of strains having a trait available


### Reading the traits tables

# Ec+Kpn strains carrying pOXA-48
traits_caps_pOXA48_refseq <- read.csv("traits_caps_pOXA-48_refseq_all.tsv", header=TRUE, sep="\t", row.names=1)
# Ec+Kpn strains carrying plasmids
traits_caps_plasmids_refseq <- read.csv("traits_caps_plasmids_refseq_all.tsv", header=TRUE, sep="\t", row.names=1)
# Ec strains carrying plasmids
traits_caps_plasmids_refseq_ec <- read.csv("traits_caps_plasmids_refseq_Ec.tsv", header=TRUE, sep="\t", row.names=1)


### Phylogenetic logistic regression

# Ec+Kpn strains carrying pOXA-48
model_caps_pOXA48_refseq <- phyloglm(Dependent~Independent, phy=tree_refseq, data=traits_caps_pOXA48_refseq, boot=100)
summary(model_caps_pOXA48_refseq)
# Ec+Kpn strains carrying plasmids
model_caps_plasmids_refseq <- phyloglm(Dependent~Independent, phy=tree_refseq, data=traits_caps_plasmids_refseq, boot=100)
summary(model_caps_plasmids_refseq)
# Ec strains carrying plasmids
model_caps_plasmids_refseq_ec <- phyloglm(Dependent~Independent, phy=tree_refseq_ec, data=traits_caps_plasmids_refseq_ec, boot=100)
summary(model_caps_plasmids_refseq_ec)

