# suppressPackageStartupMessages({
#     library(NMF)
#     library(scran)
#     library(Seurat)
#     library(SeuratData)
#     library(SingleCellExperiment)
# })
# 
# sce0 <- as.SingleCellExperiment(readRDS("~/packages/SPOTlight/sample_data/allen_cortex_dwn.rds"))
# spe0 <- as.SingleCellExperiment(LoadData("stxBrain", type = "anterior1"))
# 
# colLabels(sce0) <- sce0$subclass
# mgs0 <- scoreMarkers(sce0)
# mgs0 <- lapply(seq_along(mgs0), \(i) {
#     df <- mgs0[[i]]
#     df$gene <- rownames(df)
#     df$weight <- df$mean.AUC
#     df$cluster <- names(mgs0)[i]
#     return(df)
# })
# mgs0 <- do.call(rbind, mgs0)

mgs <- mgs0
n_top <- NULL
gene_id <- "gene"
group_id <- "cluster"
weight_id <- "mean.AUC"

sce <- sce0
spe <- sce0
groups <-  colLabels(sce, onAbsence = "error")

n_cells <- 50
n_genes <- 1e3
verbose <- TRUE

scale <- TRUE
model <- "ns"
