set.seed(321)
# mock up some single-cell, mixture & marker data
sce <- mockSC()
spe <- mockSP(sce)
mgs <- getMGS(sce)

# .scale_uv ----
test_that(".scale_uv()", {
    x <- counts(sce)
    y <- .scale_uv(x)
    expect_is(y, "matrix")
    expect_identical(dim(y), dim(x))
    expect_identical(dimnames(y), dimnames(x))
    expect_true(all(abs(1 - sparseMatrixStats::rowVars(y)) < 1e-12))
})

# default parameters
defs <- list(
    gene_id = "gene",
    group_id = "type",
    weight_id = "weight",
    verbose = FALSE
)

# NMF ----
test_that("NMF", {
    x <- counts(sce)
    y <- counts(spe)
    groups <- sce$type
    group_ids <- paste0("topic_", seq_len(length(levels(sce$type))))
    n_groups <- length(group_ids)
    
    # + trainNMF ----
    # undetected genes should be filtered out
    # and pass silently (i.e., without error)
    i <- sample(rownames(x), 5)
    j <- sample(rownames(y), 5)
    x. <- x
    x.[i, ] <- 0
    y. <- y
    y.[j, ] <- 0
    args <- c(defs, list(x., y., groups, mgs))
    fit <- expect_silent(do.call(trainNMF, args))
    mod <- fit[["mod"]]
    expect_is(mod, "NMF")
    expect_true(!all(c(i, j) %in% rownames(mod)))
    # Only marker genes should be present - we don't use hvg here
    expect_true(all(rownames(mod) %in% mgs$gene))
    # valid call should give an object of class 'NMF'
    # and dimension (#genes) x (#cells) x (#groups)
    args <- c(defs, list(x, y, groups, mgs))
    fit <- expect_silent(do.call(trainNMF, args))
    mod <- fit[["mod"]]
    expect_is(mod, "NMF")
    # Remove genes since these can change depending on 
    # filtering, mgs, hvg, all 0...
    expect_identical(
        dimnames(mod)[2:3],
        c(dimnames(x)[2], list(group_ids)))
    
    # + .topic_profiles ----
    # should give a square numeric matrix
    # of dimension (#groups) x (#groups)
    ref <- .topic_profiles(mod, groups)
    expect_is(ref, "matrix")
    expect_equal(dim(ref), rep(n_groups, 2))
    expect_identical(rownames(ref), paste0("type", 1:nrow(ref)))
    expect_identical(colnames(ref), group_ids)
    
    # + .pred_prop ----
    fqs <- .pred_prop(x, mod)
    expect_is(fqs, "matrix")
    expect_true(is.numeric(x))
    expect_true(all(fqs >= 0))
    expect_equal(dim(fqs), c(n_groups, ncol(x)))
    expect_identical(dimnames(fqs), list(group_ids, colnames(x)))
    
    # + runDeconvolution ----
    # should give a numeric matrix
    # of dimension (#groups) x (#spots)
    # with proportions (i.e., values in [0, 1])
    res <- runDeconvolution(y, mod, ref)
    mat <- res[[1]]
    err <- res[[2]]
    expect_is(mat, "matrix")
    expect_true(is.numeric(err))
    expect_true(is.numeric(x))
    expect_true(all(mat >= 0 & mat <= 1))
    expect_true(all(rowSums(mat) - 1 < 1e-12))
    expect_identical(dimnames(mat), list(colnames(y), group_ids))
    expect_identical(rownames(mat), names(err))
    # actually check the estimates are legit
    # (MSE < 0.1 compared to simulated truth)
    sim <- metadata(spe)[[1]]
    mse <- mean((mat - sim)^2)
    expect_true(mse < 0.1)
})
