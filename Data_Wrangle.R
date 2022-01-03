rm(list = ls())
library(Seurat)
library(RCurl)
library(dplyr)
library(fgsea)
library(scater)
source("utils.R")
options(stringsAsFactors = F)
# Load Data ----------------------------------------------------------------------------------------------------------------------

counts <- read.table("OMIX270-20-01.csv", sep = ",", header = 1, row.names = 1)
cell_metadata <- read.table("OMIX270-20-02.csv", sep = ",", header = 1)
patient_metadata <- read.table("Patient_MetaData.csv", sep = ",", header=1)

# Convert raw counts to TPM ------------------------------------------------------------------------------------------------------

gene.sizes <- read.delim("gene_length.txt", col.names= c('Gene', 'Length'))

# match the symbols with the count matrix and TPM normalize --------------------------------------------------------------------------------------
gene.sizes <- gene.sizes[gene.sizes$Gene %in% rownames(counts),]
counts2 <- counts[rownames(counts) %in% gene.sizes$Gene,]

gene.sizes  <- gene.sizes %>%
  dplyr::slice(match(rownames(counts2), Gene))

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

counts_tpm <- as.data.frame(tpm3(counts2, gene.sizes$Length))
all_data <- counts_tpm

# Get col_data for scater object -------------------------------------------------------------------------------------------------------

final_metadata <- merge(cell_metadata, patient_metadata, by.x = 'patient', by.y = 'Patients')
rownames(final_metadata) <- final_metadata$X
col_data <- final_metadata

# mark the metabolic genes --------------------------------------------------------------------------------------------------------------------------------------------------------------------
dim(all_data)
all_data <- all_data[,colnames(all_data) %in% rownames(col_data),]
dim(all_data)
all_data <- data.matrix(all_data)
pathways <- gmtPathways("KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE

# Build scater object --------------------------------------------------------------------------------------------------------------------------------------------------------------------

sce <- SingleCellExperiment(
  assays = list(tpm=all_data,exprs=all_data),
  colData = col_data,
  rowData = row_data
)

# Save final scater object --------------------------------------------------------------------------------
saveRDS(sce,"selected_sce.rds")

