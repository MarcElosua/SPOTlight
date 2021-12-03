# mock up some single-cell, mixture & marker data
sce <- .mock_sc(ng = 200, nc = 10, nt = 3)
spe <- .mock_sp(sce)
mgs <- .get_mgs(sce)

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
  
  mtrx <- res[[1]]
  mod <- res[[2]]
  expect_is(res, "list")
  expect_is(mtrx, "matrix")
  expect_is(mod, "NMFfit")
  expect_identical(ncol(mtrx), as.integer(length(unique(sce$type)) + 1))
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
  
  mtrx <- res[[1]]
  mod <- res[[2]]
  expect_is(res, "list")
  expect_is(mtrx, "matrix")
  expect_is(mod, "NMFfit")
  expect_identical(ncol(mtrx), as.integer(length(unique(sce$type)) + 1))
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
  
  mtrx <- res[[1]]
  mod <- res[[2]]
  expect_is(res, "list")
  expect_is(mtrx, "matrix")
  expect_is(mod, "NMFfit")
  expect_identical(ncol(mtrx), as.integer(length(unique(sce$type)) + 1))
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
  
  mtrx <- res[[1]]
  mod <- res[[2]]
  expect_is(res, "list")
  expect_is(mtrx, "matrix")
  expect_is(mod, "NMFfit")
  expect_identical(ncol(mtrx), as.integer(length(unique(sce$type)) + 1))
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
  
  mtrx <- res[[1]]
  mod <- res[[2]]
  expect_is(res, "list")
  expect_is(mtrx, "matrix")
  expect_is(mod, "NMFfit")
  expect_identical(ncol(mtrx), as.integer(length(unique(sce$type)) + 1))
})