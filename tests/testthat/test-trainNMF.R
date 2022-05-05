set.seed(321)
# mock up some single-cell, mixture & marker data
sce <- mockSC(ng = 200, nc = 10, nt = 3)
spe <- mockSP(sce)
mgs <- getMGS(sce)

# Create dummy Seurat object
sec <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts(sce)))
sep <- SeuratObject::CreateSeuratObject(counts = counts(spe))

# Function to run the checks
# .checks_nmf <- function(res, sce) {
#     mod <- res[[1]]
#     mtr <- res[[2]]
#     expect_is(res, "list")
#     expect_is(mtr, "matrix")
#     expect_is(mod, "list")
#     expect_identical(ncol(mtr), length(unique(sce$type)))
#     expect_identical(sort(rownames(mtr)), sort(unique(as.character(sce$type))))
#     expect_identical(nrow(mtr), ncol(mod$w))
#     expect_identical(nrow(mtr), nrow(mod$h))
# }

# .checks_rcpp <- function(res, sce) {
#     mod <- res[[1]]
#     mtr <- res[[2]]
#     expect_is(res, "list")
#     expect_is(mtr, "matrix")
#     expect_is(mod, "list")
#     expect_identical(ncol(mtr), length(unique(sce$type)))
#     expect_identical(nrow(mtr), ncol(mod$w))
#     expect_identical(nrow(mtr), nrow(mod$h))
# }

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


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check RCPP trainNMF x, y inputs  -------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# trainNMF with SCE ----
test_that("rcpp trainNMF x SCE", {
    res <- trainNMF(
        x = sce,
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with SPE ----
test_that("rcpp trainNMF x SPE", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = spe,
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with SPE ----
test_that("rcpp trainNMF x SEC", {
    res <- trainNMF(
        x = sec,
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with SEP ----
test_that("rcpp trainNMF x SEP", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = sep,
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})


# trainNMF with sparse matrix sc ----
test_that("rcpp trainNMF x dgCMatrix SC", {
    res <- trainNMF(
        x = Matrix::Matrix(counts(sce), sparse = TRUE),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    .checks(res, sce)
})

# trainNMF with sparse matrix sp ----
test_that("rcpp trainNMF x dgCMatrix SP", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = Matrix::Matrix(counts(spe), sparse = TRUE),
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with sparse matrix sc ----
test_that("rcpp trainNMF x DelayedMatrix SC", {
    res <- trainNMF(
        x = DelayedArray::DelayedArray(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    .checks(res, sce)
})

# trainNMF with sparse matrix sp ----
test_that("rcpp trainNMF x DelayedMatrix SP", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = DelayedArray::DelayedArray(counts(sce)),
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with matrices in both ----
test_that("rcpp trainNMF x matrices", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "RcppML",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with matrices in both and HVG----
test_that("rcpp trainNMF x hvg", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "RcppML",
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
# ----  Check NMF::nmf trainNMF x, y inputs  -------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# trainNMF with SCE ----
test_that("NMF trainNMF x SCE", {
    res <- trainNMF(
        x = sce,
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with SPE ----
test_that("NMF trainNMF x SPE", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = spe,
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with SPE ----
test_that("NMF trainNMF x SEC", {
    res <- trainNMF(
        x = sec,
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with SEP ----
test_that("NMF trainNMF x SEP", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = sep,
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})


# trainNMF with sparse matrix sc ----
test_that("NMF trainNMF x dgCMatrix SC", {
    res <- trainNMF(
        x = Matrix::Matrix(counts(sce), sparse = TRUE),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    .checks(res, sce)
})

# trainNMF with sparse matrix sp ----
test_that("NMF trainNMF x dgCMatrix SP", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = Matrix::Matrix(counts(spe), sparse = TRUE),
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with sparse matrix sc ----
test_that("NMF trainNMF x DelayedMatrix SC", {
    res <- trainNMF(
        x = DelayedArray::DelayedArray(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    .checks(res, sce)
})

# trainNMF with sparse matrix sp ----
test_that("NMF trainNMF x DelayedMatrix SP", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = DelayedArray::DelayedArray(counts(sce)),
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with matrices in both ----
test_that("NMF trainNMF x matrices", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene"
    )
    
    .checks(res, sce)
})

# trainNMF with matrices in both and HVG----
test_that("NMF trainNMF x hvg", {
    res <- trainNMF(
        x = as.matrix(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        pnmf = "NMF",
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene",
        hvg = row.names(sce)[seq_len(50)]
    )
    
    .checks(res, sce)
})

