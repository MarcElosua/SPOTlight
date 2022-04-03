set.seed(321)
# mock up some single-cell, mixture & marker data
sce <- mockSC(ng = 200, nc = 10, nt = 3)
spe <- mockSP(sce)
mgs <- getMGS(sce)
# Create SpatialExperiment
spe1 <- SpatialExperiment::SpatialExperiment(
    assay = list(counts = counts(spe)), colData = colData(spe))

# Create dummy Seurat object
sec <- SeuratObject::CreateSeuratObject(counts = counts(sce))
sep <- SeuratObject::CreateSeuratObject(counts = counts(spe))

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
# runDeconvolution with SCE ----
test_that("runDeconvolution x SCE", {
    decon <- runDeconvolution(
        x = spe,
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

test_that("runDeconvolution x SPE", {
    decon <- runDeconvolution(
        x = spe1,
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

# runDeconvolution with Seurat ----
test_that("runDeconvolution x SEP", {
    decon <- runDeconvolution(
        x = sep,
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

# runDeconvolution with sparse matrix sp ----
test_that("runDeconvolution x dgCMatrix SP", {
    decon <- runDeconvolution(
        x = Matrix::Matrix(counts(spe), sparse = TRUE),
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

# runDeconvolution with sparse matrix sp ----
test_that("runDeconvolution x DelayedMatrix SP", {
    decon <- runDeconvolution(
        x = DelayedArray::DelayedArray(counts(sce)),
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

# runDeconvolution with matrices in both ----
test_that("runDeconvolution x matrices", {
    decon <- runDeconvolution(
        x = as.matrix(counts(spe)),
        mod = res[["mod"]],
        ref = res[["topic"]]
    )
    
    .checks(decon, sce)
})

