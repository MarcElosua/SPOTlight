set.seed(321)
# mock up some single-cell, mixture & marker data
sce <- .mock_sc(ng = 200, nc = 10, nt = 3)
spe <- .mock_sp(sce)
mgs <- .get_mgs(sce)

# Function to run the checks
.checks <- function(res, sce) {
    mtr <- res[[1]]
    rss <- res[[2]]
    mod <- res[[3]]
    expect_is(res, "list")
    expect_is(mtr, "matrix")
    expect_is(rss, "numeric")
    expect_is(mod, "NMFfit")
    expect_identical(ncol(mtr), length(unique(sce$type)))
    expect_identical(nrow(mtr), length(rss))
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----  Check SPOTlight x, y inputs  -------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# .SPOTlight with SCE ----
test_that("SPOTlight x SCE", {
  
    res <- SPOTlight(
        x = sce,
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene")
  
    .checks(res, sce)
})

# .SPOTlight with SPE ----
test_that("SPOTlight x SPE", {
  
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = spe,
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene")
  
    .checks(res, sce)
    })


# .SPOTlight with sparse matrix sc ----
test_that("SPOTlight x sparse SC", {
  
    res <- SPOTlight(
        x = counts(sce),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene")
    .checks(res, sce)
    })

# .SPOTlight with sparse matrix sp ----
test_that("SPOTlight x sparse SP", {
  
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = counts(spe),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene")
  
    .checks(res, sce)  
})

# .SPOTlight with matrices in both ----
test_that("SPOTlight x sparse SP", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene")
  
    .checks(res, sce)
})

# .SPOTlight with matrices in both and HVG----
test_that("SPOTlight x sparse SP", {
    res <- SPOTlight(
        x = as.matrix(counts(sce)),
        y = as.matrix(counts(spe)),
        groups = sce$type,
        mgs = mgs,
        weight_id = "weight",
        group_id = "type",
        gene_id = "gene",
        hvg = row.names(sce)[1:50])
    
    .checks(res, sce)
})