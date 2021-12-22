set.seed(321)
x <- replicate(m <- 25, runif(10, 0, 1))
# Add an anticorrelated column
x[, 24] <- seq(0, 1, length.out = 10)
x[, 25] <- seq(1, 0, length.out = 10)
rownames(x) <- paste0("spot", 1:nrow(x))
colnames(x) <- paste0("type", 1:ncol(x))

.checks <- function(p) {
    expect_is(p, "ggplot")
    expect_true(all(p$data$p >= 0))
    expect_true(all(p$data$p <= 1))
    expect_true(is.numeric(p$data$value))
    expect_true(max(p$data$value) == 1)
    expect_true(nrow(p$data) == m * m)
}

# plotCorrelationMatrix basic ----
test_that("plotCorrelationMatrix basic", {
    # The most basic example
    p <- plotCorrelationMatrix(x = x)
    .checks(p)
})

# plotCorrelationMatrix() spearman correlation ----
test_that("plotCorrelationMatrix() spearman", {
    # The most basic example
    p <- plotCorrelationMatrix(
        x = x,
        cor.method = "kendall"
    )
    .checks(p)
})

# plotCorrelationMatrix() ----
test_that("plotCorrelationMatrix() insig", {
    # The most basic example
    p <- plotCorrelationMatrix(
        x = x,
        insig = "pch"
    )

    .checks(p)
    # This adds an extra layer with the X on top of the insig
    expect_true(length(p$layers) == 2)
})

# plotCorrelationMatrix() colors ----
test_that("plotCorrelationMatrix() colors", {
    # The most basic example
    p <- plotCorrelationMatrix(
        x = x,
        colors = c("#FF00FF", "#FFFFFF", "#000000")
    )

    .checks(p)
    g <- ggplot_build(p)

    # max color
    i <- which(p$data$value == max(p$data$value))[[1]]
    expect_identical(g$data[[1]][i, ][, "fill"], "#000000")
    # 0 color
    j <- which(p$data$value == 0)[[1]]
    expect_identical(g$data[[1]][j, ][, "fill"], "#FFFFFF")

    # min color
    k <- which(p$data$value == min(p$data$value))[[1]]
    expect_identical(g$data[[1]][k, ][, "fill"], "#FF00FF")
})

# plotCorrelationMatrix() hc.order ----
test_that("plotCorrelationMatrix() hc.order", {
    # The most basic example
    p <- plotCorrelationMatrix(
        x = x,
        hc.order = FALSE
    )

    .checks(p)
    # Make sure the order is no changed
    expect_equal(as.character(p$data$Var1[1:ncol(x)]), colnames(x))
})

# plotCorrelationMatrix() p.mat ----
test_that("plotCorrelationMatrix() p.mat", {
    # The most basic example
    p <- plotCorrelationMatrix(
        x = x,
        p.mat = FALSE
    )

    # Make sure the p value is not computed
    expect_is(p, "ggplot")
    expect_true(is.numeric(p$data$value))
    expect_true(max(p$data$value) == 1)
    expect_true(nrow(p$data) == m * m)
    expect_true(all(is.na(p$data$pvalue)))
    expect_true(all(is.na(p$data$signif)))
})
