# helper to record base R plot
. <- \(.) {
    pdf(NULL)
    dev.control(displaylist = "enable")
    set.seed(1); .
    . <- recordPlot()
    invisible(dev.off())
    return(.)
}

# mock up some data
set.seed(321)
x <- replicate(m <- 10, rnorm(n <- 100, runif(1, -1, 1)))
y <- replicate(m <- 10, rnorm(n <- 100, runif(1, -1, 1)))

test_that("plotInteractions(), which = 'heatmap'", {
    p <- plotInteractions(x, "heatmap")
    expect_is(p, "ggplot")
    expect_true(all(p$data$p >= 0))
    expect_true(all(p$data$p <= 1))
    expect_true(is.integer(p$data$n))
    expect_true(nrow(p$data) == m*(m-1)/2)
})

test_that("plotInteractions(), which = 'network'", {
    # TODO
    expect_true(TRUE)
})
