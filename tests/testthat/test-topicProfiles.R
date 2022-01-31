set.seed(123)
x <- mockSC(nc = 50, nt = 3)
y <- mockSP(x)
z <- getMGS(x)
res <- SPOTlight(x, y,
    groups = x$type,
    mgs = z,
    group_id = "type",
    verbose = FALSE)

test_that("plotTopicProfiles common", {
    p <- plotTopicProfiles(res[[3]], x$type, facet = FALSE)
    expect_is(p, "ggplot")
    expect_equal(nrow(p$data), 9)
    expect_equal(ncol(p$data), 3)
    g <- ggplot_build(p)
    expect_equal(unique(g$data[[1]]$colour), c("#0000FF", "#D3D3D3"))
})

test_that("plotTopicProfiles facet", {
    p <- plotTopicProfiles(res[[3]], x$type, facet = TRUE)
    expect_is(p, "ggplot")
    expect_equal(nrow(p$data), 150)
    expect_equal(ncol(p$data), 4)
    g <- ggplot_build(p)
    expect_equal(unique(g$data[[1]]$colour), c("#0000FF", "#D3D3D3"))
})
