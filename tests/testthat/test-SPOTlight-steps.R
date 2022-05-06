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

###############################
#### Run SPOTlight wrapper ####
###############################
set.seed(687)
res1 <- SPOTlight(
    x = sce,
    y = as.matrix(counts(spe)),
    groups = sce$type,
    mgs = mgs,
    weight_id = "weight",
    group_id = "type",
    gene_id = "gene"
)

################################
#### Run SPOTlight by steps ####
################################
set.seed(687)
# Train NMF
mod_ls <- trainNMF(
    x = as.matrix(counts(sce)),
    y = as.matrix(counts(spe)),
    groups = sce$type,
    mgs = mgs,
    weight_id = "weight",
    group_id = "type",
    gene_id = "gene"
)

res2 <- runDeconvolution(
    x = spe,
    mod = mod_ls[["mod"]],
    ref = mod_ls[["topic"]]
)

# NMF ----
test_that("SPOTlight vs SPOTlight-steps", {

    # basis and coef should be the same between SPOTlight and SPOTlight-steps
    expect_true(all(basis(res1[["NMF"]]) == basis(mod_ls[["mod"]])))
    expect_true(all(coef(res1[["NMF"]]) == coef(res2[["NMF"]])))
    
    # Deconvolution results are the same
    expect_true(all(res1[["mat"]] == res2[["mat"]]))
    
    # actually check the estimates are legit
    # (MSE < 0.1 compared to simulated truth)
    sim <- S4Vectors::metadata(spe)[[1]]
    mse <- mean((res2[["mat"]] - sim)^2)
    expect_true(mse < 0.1)
    
})

