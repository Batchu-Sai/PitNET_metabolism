rm(list = ls())
library(scater)
library(reshape2)

# Loading the tumor data ------------------------------------------------------------------
selected_sce <- readRDS("selected_sce.rds")
selected_metabolic_sce <- selected_sce[rowData(selected_sce)$metabolic,]

# T-SNE ===========================================================================
library("Rtsne")
set.seed(12345)

tsne_metabolic <- Rtsne(t(assay(selected_metabolic_sce,"exprs")),initial_dims=20,theta=0.0,perplexity = 30)
tsne_metabolic_out <- data.frame(x=tsne_metabolic$Y[,1],y=tsne_metabolic$Y[,2],group = colData(selected_metabolic_sce)$X2017_WHO, lineage = colData(selected_metabolic_sce)$lineage)

# Visualize
ggplot(tsne_metabolic_out) + geom_point(aes(x, y, colour = group), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Using metabolic genes across 2017 WHO")

ggplot(tsne_metabolic_out) + geom_point(aes(x, y, colour = lineage), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Using metabolic genes across lineages")

