set.seed(321)
# mock up some single-cell, mixture & marker data
sce <- mockSC(ng = 200, nc = 10, nt = 3)
spe <- mockSP(sce)
mgs <- getMGS(sce)

.checks <- function(res, sce) {
    mod <- res[[1]]
    mtr <- res[[2]]
    expect_is(res, "list")
    expect_is(mtr, "matrix")
    expect_is(mod, "list")
    expect_identical(ncol(mtr), length(unique(sce$type)))
    expect_identical(nrow(mtr), ncol(mod$w))
    expect_identical(nrow(mtr), nrow(mod$h))
    }

# Unit test to verify that the topic with max weight in each cell aligns with type
.check_topic_alignment <- function(mod_h, topic) {
    cell_pos <- c(1, 11, 21)
    type_names <- rownames(topic)[1:3]   # types 1 to 3
    
    for (i in seq_along(cell_pos)) {
        cell <- cell_pos[i]
        type <- type_names[i]
        expect_equal(which.max(mod_h[, cell]), which.max(topic[type, ]))
    }
}
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check RCPP trainNMF x, y inputs  -------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# trainNMF with SCE ----
test_that("rcpp trainNMF x SCE", {
    set.seed(321)
    res <- trainNMF(
        x = sce,
        y = rownames(spe),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
    .check_topic_alignment(res$mod$h, res$topic)
})


# trainNMF with sparse matrix sc ----
test_that("rcpp trainNMF x dgCMatrix SC", {
    set.seed(321)
    res <- trainNMF(
        x = Matrix::Matrix(counts(sce), sparse = TRUE),
        y = rownames(spe),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    .checks(res, sce)
    .check_topic_alignment(res$mod$h, res$topic)
})

# trainNMF with sparse matrix sc ----
test_that("rcpp trainNMF x DelayedMatrix SC", {
    set.seed(321)
    res <- trainNMF(
        x = DelayedArray::DelayedArray(counts(sce)),
        y = rownames(spe),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    .checks(res, sce)
    .check_topic_alignment(res$mod$h, res$topic)
})

# trainNMF with matrices in both ----
test_that("rcpp trainNMF x matrices", {
    set.seed(321)
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = rownames(spe),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
    .check_topic_alignment(res$mod$h, res$topic)
})

# trainNMF with matrices in both and HVG----
test_that("rcpp trainNMF x hvg", {
    set.seed(321)
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = rownames(spe),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene",
        hvg = row.names(sce)[seq_len(50)]
    )
    
    .checks(res, sce)
    .check_topic_alignment(res$mod$h, res$topic)
})


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check NMF::nmf trainNMF x, y inputs  -------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# trainNMF with SCE ----
# test_that("NMF trainNMF x SCE", {
#     res <- trainNMF(
#         x = sce,
#         y = rownames(spe),
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
# # trainNMF with SPE ----
# test_that("NMF trainNMF x SEC", {
#     res <- trainNMF(
#         x = sec,
#         y = rownames(spe),
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
# # trainNMF with sparse matrix sc ----
# test_that("NMF trainNMF x dgCMatrix SC", {
#     res <- trainNMF(
#         x = Matrix::Matrix(counts(sce), sparse = TRUE),
#         y = rownames(spe),
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
# # trainNMF with sparse matrix sc ----
# test_that("NMF trainNMF x DelayedMatrix SC", {
#     res <- trainNMF(
#         x = DelayedArray::DelayedArray(counts(sce)),
#         y = rownames(spe),
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
# # trainNMF with matrices in both ----
# test_that("NMF trainNMF x matrices", {
#     res <- trainNMF(
#         x = as.matrix(counts(sce)),
#         y = rownames(spe),
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
# # trainNMF with matrices in both and HVG----
# test_that("NMF trainNMF x hvg", {
#     res <- trainNMF(
#         x = as.matrix(counts(sce)),
#         y = rownames(spe),
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
