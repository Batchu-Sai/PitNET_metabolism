rm(list = ls())
library(scater)
library(stringr)
library(pheatmap)
library(gtools)
library(scran)
source("utils.R")
options(stringsAsFactors=FALSE)
source("runGSEA_preRank.R")
pathway_file <- "KEGG_metabolism.gmt"

# Loading the data and rename tumors to prevent double names or spaces

selected_tumor_sce <- readRDS("selected_sce.rds")
colData(selected_tumor_sce)[colData(selected_tumor_sce)$X2017_WHO == 'Gonadotroph adenoma',]$X2017_WHO <- "Gonado"
colData(selected_tumor_sce)[colData(selected_tumor_sce)$X2017_WHO == 'Lactotroph adenomas',]$X2017_WHO <- "Lacto"
colData(selected_tumor_sce)[colData(selected_tumor_sce)$X2017_WHO == 'ACTH adenomas',]$X2017_WHO <- "Cortico"
colData(selected_tumor_sce)[colData(selected_tumor_sce)$X2017_WHO == 'Somatotroph adenomas',]$X2017_WHO <- "Somato"
colData(selected_tumor_sce)[colData(selected_tumor_sce)$X2017_WHO == 'PIT1+ PPA',]$X2017_WHO <- "PIT1_PPA"
selected_tumor_metabolic_sce <- selected_tumor_sce[rowData(selected_tumor_sce)$metabolic,]
#=========================================================================
tumors <- unique(selected_tumor_sce$X2017_WHO)

# Tumor cells
enrich_data_df <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)

pc_plotdata <- data.frame(x=numeric(),y=numeric(),
                          sel=character(),tumor=character())

for (t in tumors){
  each_metabolic_sce <- selected_tumor_metabolic_sce[,selected_tumor_metabolic_sce$X2017_WHO==t]
  each_metabolic_tpm <- assay(each_metabolic_sce,"exprs")
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]
  x <- as.matrix(each_metabolic_tpm)
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  
  ###select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>=0.8)[1]
  
  ###plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             tumor=rep(t,length(percentVar)))
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)
  
  ####
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  
  runGSEA_preRank(pre_rank_matrix,pathway_file,t)
  #get the result
  result_dir <- list.files(path="preRankResults",pattern = paste0("^",t,".GseaPreranked(.*)"),full.names=T)
  result_file <- list.files(path=result_dir,pattern="gsea_report_for_na_pos_(.*).xls",full.names=T)
  gsea_result <- read.table(result_file,header = T,sep="\t",row.names=1)
  gsea_pathways <- str_to_title(rownames(gsea_result))
  gsea_pathways <- str_replace(gsea_pathways,"Tca","TCA")
  gsea_pathways <- str_replace(gsea_pathways,"Gpi","GPI")
  enrich_data_df <- rbind(enrich_data_df,data.frame(x=t,y=gsea_pathways,NES=gsea_result$NES,PVAL=gsea_result$NOM.p.val))
}

write.csv(enrich_data_df, quote = F, "enriched_data_group.csv", row.names = F)

##plot variance
pc_plotdata$tumor <- factor(pc_plotdata$tumor,levels=mixedsort(tumors))
p <- ggplot(pc_plotdata) + geom_point(aes(x,y,colour=factor(sel)),size=0.5) +
  scale_color_manual(values=c("gray","#ff4000")) +
  facet_wrap(~tumor,scales="free",ncol = 4) + theme_bw() + 
  labs(x="Principal components", y="Explained variance (%)") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor= element_blank(),
        axis.line=element_line(size=0.2,colour="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(colour="black", size = 6),
        strip.background = element_rect(fill="white",size=0.2,colour = NULL),
        strip.text=element_text(size=6))

ggsave(file.path(outDir,"malignant_PC_variance_plot.pdf"),p,device="pdf",useDingbats=FALSE)
unlink("preRankResults",recursive=T)
unlink("prerank.rnk")
date_string <- Sys.Date()
date_split <- strsplit(as.character(date_string),"-")[[1]]
unlink(paste0(tolower(month.abb[as.numeric(date_split[2])]),date_split[3]),recursive=T)




# For Lineage ===========================================================================================================
lineage <- unique(selected_tumor_sce$lineage)
enrich_data_df <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)
pc_plotdata <- data.frame(x=numeric(),y=numeric(),
                          sel=character(),tumor=character())

for (t in lineage){
  each_metabolic_sce <- selected_tumor_metabolic_sce[,selected_tumor_metabolic_sce$lineage==t]
  each_metabolic_tpm <- assay(each_metabolic_sce,"exprs")
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]
  x <- as.matrix(each_metabolic_tpm)
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  
  ###select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>=0.8)[1]
  
  ###plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             tumor=rep(t,length(percentVar)))
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)
  
  ####
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  
  runGSEA_preRank(pre_rank_matrix,pathway_file,t)
  #get the result
  result_dir <- list.files(path="preRankResults",pattern = paste0("^",t,".GseaPreranked(.*)"),full.names=T)
  result_file <- list.files(path=result_dir,pattern="gsea_report_for_na_pos_(.*).xls",full.names=T)
  gsea_result <- read.table(result_file,header = T,sep="\t",row.names=1)
  gsea_pathways <- str_to_title(rownames(gsea_result))
  gsea_pathways <- str_replace(gsea_pathways,"Tca","TCA")
  gsea_pathways <- str_replace(gsea_pathways,"Gpi","GPI")
  enrich_data_df <- rbind(enrich_data_df,data.frame(x=t,y=gsea_pathways,NES=gsea_result$NES,PVAL=gsea_result$NOM.p.val))
}
write.csv(enrich_data_df, quote = F, "enriched_data_lineage.csv", row.names = F)

##plot variance
pc_plotdata$tumor <- factor(pc_plotdata$tumor,levels=mixedsort(lineage))
p <- ggplot(pc_plotdata) + geom_point(aes(x,y,colour=factor(sel)),size=0.5) +
  scale_color_manual(values=c("gray","#ff4000")) +
  facet_wrap(~tumor,scales="free",ncol = 4) + theme_bw() + 
  labs(x="Principal components", y="Explained variance (%)") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor= element_blank(),
        axis.line=element_line(size=0.2,colour="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(colour="black", size = 6),
        strip.background = element_rect(fill="white",size=0.2,colour = NULL),
        strip.text=element_text(size=6))

ggsave(file.path(outDir,"malignant_PC_variance_plot.pdf"),p,device="pdf",useDingbats=FALSE)
unlink("preRankResults",recursive=T)
unlink("prerank.rnk")
date_string <- Sys.Date()
date_split <- strsplit(as.character(date_string),"-")[[1]]
unlink(paste0(tolower(month.abb[as.numeric(date_split[2])]),date_split[3]),recursive=T)