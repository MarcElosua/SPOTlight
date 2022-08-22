set.seed(321)
# mock up some single-cell, mixture & marker data
sce <- mockSC(ng = 200, nc = 10, nt = 3)
spe <- mockSP(sce)
mgs <- getMGS(sce)
# Create SpatialExperiment
spe1 <- SpatialExperiment::SpatialExperiment(
    assay = list(counts = counts(spe)),
    colData = SummarizedExperiment::colData(spe))

# Create dummy Seurat object
sec <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts(sce)))
sep <- SeuratObject::CreateSeuratObject(counts = counts(spe))

# Function to run the checks
.checks <- function(res, sce) {
    mtr <- res[[1]]
    rss <- res[[2]]
    mod <- res[[3]]
    expect_is(res, "list")
    expect_is(mtr, "matrix")
    expect_is(rss, "numeric")
    expect_is(mod, "list")
    expect_identical(ncol(mtr), length(unique(sce$type)))
    expect_identical(sort(colnames(mtr)), sort(unique(as.character(sce$type))))
    expect_identical(nrow(mtr), length(rss))
    expect_identical(sort(rownames(mtr)), sort(names(rss)))
}

# .checks <- function(res, sce) {
#     mod <- res[[1]]
#     mtr <- res[[2]]
#     expect_is(res, "list")
#     expect_is(mtr, "matrix")
#     expect_is(mod, "list")
#     expect_identical(ncol(mtr), length(unique(sce$type)))
#     expect_identical(nrow(mtr), ncol(mod$w))
#     expect_identical(nrow(mtr), nrow(mod$h))
# }

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check SPOTlight x, y inputs  -------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# .SPOTlight with SCE ----
test_that("SPOTlight x SCE rcpp", {
    res <- SPOTlight(
        x = sce,
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# .SPOTlight with SPE ----
test_that("SPOTlight x SCE spatial rcpp", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = spe,
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )

    .checks(res, sce)
})

# .SPOTlight with SPE ----
test_that("SPOTlight x SCE spatial rcpp", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = spe1,
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# .SPOTlight with Seurat SC ----
test_that("SPOTlight x SEC rcpp", {
    res <- SPOTlight(
        x = sec,
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# .SPOTlight with Seurat SP ----
test_that("SPOTlight x SEP rcpp", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = spe,
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# .SPOTlight with sparse matrix sc ----
test_that("SPOTlight x dgCMatrix SC rcpp", {
    res <- SPOTlight(
        x = Matrix::Matrix(counts(sce), sparse = TRUE),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    .checks(res, sce)
})

# .SPOTlight with sparse matrix sp ----
test_that("SPOTlight x dgCMatrix SP rcpp", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = Matrix::Matrix(counts(spe), sparse = TRUE),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )

    .checks(res, sce)
})

# .SPOTlight with sparse matrix sc ----
test_that("SPOTlight x DelayedMatrix SC rcpp", {
    res <- SPOTlight(
        x = DelayedArray::DelayedArray(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    .checks(res, sce)
})

# .SPOTlight with sparse matrix sp ----
test_that("SPOTlight x DelayedMatrix SP rcpp", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = DelayedArray::DelayedArray(counts(sce)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )

    .checks(res, sce)
})

# .SPOTlight with matrices in both ----
test_that("SPOTlight x matrices rcpp", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )

    .checks(res, sce)
})

# .SPOTlight with matrices in both and HVG----
test_that("SPOTlight x hvg rcpp", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene",
        hvg = row.names(sce)[seq_len(50)]
    )

    .checks(res, sce)
})


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check SPOTlight x, y inputs with NMF  ----------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# .SPOTlight with SCE ----
# test_that("SPOTlight x SCE NMF", {
#     res <- SPOTlight(
#         x = sce,
#         y = as.matrix(counts(spe)),
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     
#     .checks(res, sce)
# })
# 
# # .SPOTlight with SPE ----
# test_that("SPOTlight x SCE spatial NMF", {
#     res <- SPOTlight(
#         x = as.matrix(counts(sce)),
#         y = spe,
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     
#     .checks(res, sce)
# })
# 
# # .SPOTlight with SPE ----
# test_that("SPOTlight x SCE spatial NMF", {
#     res <- SPOTlight(
#         x = as.matrix(counts(sce)),
#         y = spe1,
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     
#     .checks(res, sce)
# })
# 
# # .SPOTlight with Seurat SC ----
# test_that("SPOTlight x SEC NMF", {
#     res <- SPOTlight(
#         x = sec,
#         y = as.matrix(counts(spe)),
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     
#     .checks(res, sce)
# })
# 
# # .SPOTlight with Seurat SP ----
# test_that("SPOTlight x SEP NMF", {
#     res <- SPOTlight(
#         x = as.matrix(counts(sce)),
#         y = spe,
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     
#     .checks(res, sce)
# })
# 
# # .SPOTlight with sparse matrix sc ----
# test_that("SPOTlight x dgCMatrix SC NMF", {
#     res <- SPOTlight(
#         x = Matrix::Matrix(counts(sce), sparse = TRUE),
#         y = as.matrix(counts(spe)),
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     .checks(res, sce)
# })
# 
# # .SPOTlight with sparse matrix sp ----
# test_that("SPOTlight x dgCMatrix SP NMF", {
#     res <- SPOTlight(
#         x = as.matrix(counts(sce)),
#         y = Matrix::Matrix(counts(spe), sparse = TRUE),
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     
#     .checks(res, sce)
# })
# 
# # .SPOTlight with sparse matrix sc ----
# test_that("SPOTlight x DelayedMatrix SC NMF", {
#     res <- SPOTlight(
#         x = DelayedArray::DelayedArray(counts(sce)),
#         y = as.matrix(counts(spe)),
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     .checks(res, sce)
# })
# 
# # .SPOTlight with sparse matrix sp ----
# test_that("SPOTlight x DelayedMatrix SP NMF", {
#     res <- SPOTlight(
#         x = as.matrix(counts(sce)),
#         y = DelayedArray::DelayedArray(counts(sce)),
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     
#     .checks(res, sce)
# })
# 
# # .SPOTlight with matrices in both ----
# test_that("SPOTlight x matrices NMF", {
#     res <- SPOTlight(
#         x = as.matrix(counts(sce)),
#         y = as.matrix(counts(spe)),
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene"
#     )
#     
#     .checks(res, sce)
# })
# 
# # .SPOTlight with matrices in both and HVG----
# test_that("SPOTlight x hvg NMF", {
#     res <- SPOTlight(
#         x = as.matrix(counts(sce)),
#         y = as.matrix(counts(spe)),
#         groups = sce$type,
#         pnmf = "NMF",
#         mgs = mgs,
#         weight_id = "weight",
#         group_id = "type",
#         gene_id = "gene",
#         hvg = row.names(sce)[seq_len(50)]
#     )
#     
#     .checks(res, sce)
# })
# 
