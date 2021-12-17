set.seed(321)
x <- replicate(m <- 10, rnorm(n <- 100, runif(1, -1, 1)))
y <- replicate(m <- 10, rnorm(n <- 100, runif(1, -1, 1)))

# plotInteractions() ----
test_that("plotInteractions() - heatmap", {
    p <- plotInteractions(x, "heatmap")
    expect_is(p, "ggplot")
    expect_true(all(p$data$p >= 0))
    expect_true(all(p$data$p <= 1))
    expect_true(is.integer(p$data$n))
    expect_true(nrow(p$data) == m*(m-1)/2)
})

# plotNetwork() ----
# record base R plot
. <- \(.) {
    pdf(NULL)
    dev.control(displaylist = "enable")
    set.seed(1); .
    . <- recordPlot()
    invisible(dev.off())
    return(.)
}

test_that("plotInteractions() - network", {
    p <- .(plotInteractions(x, "network"))
    expect_is(p, "recordedplot")
    recordedplot})
    