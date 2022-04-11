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
    p <- plotTopicProfiles(x = res[[3]], y = x$type, facet = FALSE)
    expect_is(p, "ggplot")
    expect_equal(nrow(p$data), 9)
    expect_equal(sort(unique(p$data$group)), as.character(sort(unique(x$type))))
    expect_equal(ncol(p$data), 3)
    expect_equal(
        as.character(sort(unique(p$data$topic))),
        as.character(seq_len(length(unique(x$type)))))
    g <- ggplot_build(p)
    expect_equal(unique(g$data[[1]]$colour), c("#3D2BFF", "#D3D3D3"))
})

test_that("plotTopicProfiles facet", {
    p <- plotTopicProfiles(res[[3]], x$type, facet = TRUE)
    expect_is(p, "ggplot")
    expect_equal(nrow(p$data), 151)
    expect_equal(ncol(p$data), 4)
    g <- ggplot_build(p)
    expect_equal(unique(g$data[[1]]$colour), c("#3D2BFF", "#4931FE", "#D3D3D3"))
})

