# suppressPackageStartupMessages({
#     library(Matrix)
#     library(data.table)
#     library(Seurat)
#     library(SeuratData)
#     library(dplyr)
#     library(gt)
#     library(SPOTlight)
#     library(igraph)
#     library(RColorBrewer)
# })
# 
# # load data ----
# 
# path_to_data <- system.file(package = "SPOTlight")
# cortex_sc <- readRDS(glue::glue("{path_to_data}/allen_cortex_dwn.rds"))
# anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")
# 
# # pre-processing ----
# 
# set.seed(123)
# cortex_sc <- cortex_sc %>% 
#     Seurat::SCTransform(., verbose = FALSE) %>%
#     Seurat::RunPCA(., verbose = FALSE) %>%
#     Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)
# 
# # marker genes ----
# 
# Seurat::Idents(object = cortex_sc) <- cortex_sc@meta.data$subclass
# cluster_markers_all <- 
#     Seurat::FindAllMarkers(
#         object = cortex_sc, 
#         assay = "SCT",
#         slot = "data",
#         verbose = TRUE, 
#         only.pos = TRUE)
# 
# # decomposition ----
# 
# set.seed(123)
# 
# spotlight_ls <- spotlight_deconvolution(
#     se_sc = cortex_sc,
#     counts_spatial = anterior@assays$Spatial@counts,
#     clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
#     cluster_markers = cluster_markers_all, # Dataframe with the marker genes
#     cl_n = 100, # number of cells per cell type to use
#     hvg = 3000, # Number of HVG to use
#     ntop = NULL, # How many of the marker genes to use (by default all)
#     transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
#     method = "nsNMF", # Factorization method
#     min_cont = 0) # Remove those cells contributing to a spot below a certain threshold
# 
# # assessment ----
# 
# h <- NMF::coef(nmf_mod[[1]])
# rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
# topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
#     h = h,
#     train_cell_clust = nmf_mod[[2]])
# 
# topic_profile_plts[[2]] + ggplot2::theme(
#     axis.text.x = ggplot2::element_text(angle = 90),
#     axis.text = ggplot2::element_text(size = 12))
# 
# topic_profile_plts[[1]] + theme(
#     axis.text.x = element_text(angle = 90), 
#     axis.text = element_text(size = 12))
# 
# basis_spotlight <- data.frame(NMF::basis(nmf_mod[[1]]))
# 
# colnames(basis_spotlight) <- unique(stringr::str_wrap(nmf_mod[[2]], width = 30))
# 
# basis_spotlight %>%
#     dplyr::arrange(desc(Astro)) %>%
#     round(., 5) %>%
#     DT::datatable(., filter = "top")
# 
# # visualization ----
# 
# # This is the equivalent to setting min_cont to 0.04
# decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
# decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
# decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
# rownames(decon_mtrx) <- colnames(anterior)
# 
# decon_df <- decon_mtrx %>%
#     data.frame() %>%
#     tibble::rownames_to_column("barcodes")
# 
# anterior@meta.data <- anterior@meta.data %>%
#     tibble::rownames_to_column("barcodes") %>%
#     dplyr::left_join(decon_df, by = "barcodes") %>%
#     tibble::column_to_rownames("barcodes")
# 
# Seurat::SpatialFeaturePlot(
#     object = anterior,
#     features = c("L2.3.IT", "L6b", "Meis2", "Oligo"),
#     alpha = c(0.1, 1))
# 
# cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
# 
# SPOTlight::spatial_scatterpie(
#     se_obj = anterior,
#     cell_types_all = cell_types_all,
#     img_path = "sample_data/spatial/tissue_lowres_image.png",
#     pie_scale = 0.4)
# 
# SPOTlight::spatial_scatterpie(
#     se_obj = anterior,
#     cell_types_all = cell_types_all,
#     img_path = "sample_data/spatial/tissue_lowres_image.png",
#     cell_types_interest = "L6b",
#     pie_scale = 0.8)