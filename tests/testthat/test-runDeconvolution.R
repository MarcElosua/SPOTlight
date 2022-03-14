set.seed(321)
# mock up some single-cell, mixture & marker data
sce <- mockSC(ng = 200, nc = 10, nt = 3)
spe <- mockSP(sce)
mgs <- getMGS(sce)

# Function to run the checks
.checks <- function(decon, sce) {
    mtr <- decon[[1]]
    rss <- decon[[2]]
    expect_is(decon, "list")
    expect_is(mtr, "matrix")
    expect_is(rss, "numeric")
    expect_identical(ncol(mtr), length(unique(sce$type)))
    expect_identical(nrow(mtr), length(rss))
}

# Train NMF
res <- trainNMF(
    x = as.matrix(counts(sce)),
    y = as.matrix(counts(spe)),
    groups = sce$type,
    mgs = mgs,
    weight_id = "weight",
    group_id = "type",
    gene_id = "gene"
)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check SPOTlight x, y inputs  -------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# .SPOTlight with SCE ----
test_that("SPOTlight x SPE", {
    decon <- runDeconvolution(
        x = spe,
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

# .SPOTlight with sparse matrix sp ----
test_that("SPOTlight x dgCMatrix SP", {
    decon <- runDeconvolution(
        x = Matrix::Matrix(counts(spe), sparse = TRUE),
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

# .SPOTlight with sparse matrix sp ----
test_that("SPOTlight x DelayedMatrix SP", {
    decon <- runDeconvolution(
        x = DelayedArray::DelayedArray(counts(sce)),
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

# .SPOTlight with matrices in both ----
test_that("SPOTlight x matrices", {
    decon <- runDeconvolution(
        x = as.matrix(counts(spe)),
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

