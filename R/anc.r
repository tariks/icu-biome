library(ANCOMBC)
hist(trees$Height, breaks = 10, col = "orange")
library(phyloseq)

x <- as.matrix(read.csv('./genus80_top30.csv',row.names=1))
x <- as.matrix(read.csv('../feature_tables/52N_genus_t80_selbal.csv',row.names=1))
otu <- otu_table(x,taxa_are_rows = FALSE)

meta <- as.data.frame(read.csv('../meta/meta52_current.csv',row.names=1))
samples <- sample_data(meta)
phylo <- phyloseq(otu,samples)
result <- ancombc(phylo, "M + Age", conserve=TRUE)
result$res
result$res_global
names(result)
