# helper to record base R plot
# plotNetwork() ----
# record base R plot
. <- \(.) {
    pdf(NULL)
    dev.control(displaylist = "enable")
    set.seed(1)
    .
    . <- recordPlot()
    invisible(dev.off())
    return(.)
}

# mock up some data
set.seed(321)
x <- replicate(m <- 10, rnorm(n <- 100, runif(1, -1, 1)))
# Add columns to check characteristics of interest
x[, 8] <- 1
x[, 9] <- 1
x[, 10] <- -1

test_that("plotInteractions(), which = 'heatmap', metric = default", {
    p <- plotInteractions(x, "heatmap")
    expect_is(p, "ggplot")
    expect_true(all(p$data$p >= 0))
    expect_true(all(p$data$p <= 1))
    expect_true(is.integer(p$data$n))
    expect_true(nrow(p$data) == m * (m - 1) / 2)
})

test_that("plotInteractions(), which = 'heatmap', metric = 'jaccard'", {
    p <- plotInteractions(x, "heatmap", "jaccard")
    expect_is(p, "ggplot")
    expect_true(all(p$data$p >= 0))
    expect_true(all(p$data$p <= 1))
    expect_true(is.integer(p$data$n))
    expect_true(nrow(p$data) == m * (m - 1) / 2)
})

test_that("plotInteractions(), which = 'heatmap', tunning", {
    p <- plotInteractions(x, "heatmap") +
        scale_fill_gradient(low = "#FFFF00", high = "#FF0000")

    # Same base checks
    expect_is(p, "ggplot")
    expect_true(all(p$data$p >= 0))
    expect_true(all(p$data$p <= 1))
    expect_true(is.integer(p$data$n))
    expect_true(nrow(p$data) == m * (m - 1) / 2)

    # Color checks
    g <- ggplot_build(p)
    d1 <- g$data[[1]]
    d2 <- g$data[[2]]
    # Access through tiles coordinates
    min <- d1[d1$x == 10 & d1$y == 9, "fill"]
    expect_equal(min, "#FFFF00")
    max <- d1[d1$x == 7 & d1$y == 1, "fill"]
    expect_equal(max, "#FF0000")
    na <- d2[d2$x == 2 & d2$y == 1, "fill"]
    expect_equal(na, "grey50")
})


test_that("plotInteractions(), which = 'network', metric = 'default'", {
    p <- .(plotInteractions(x, "network"))
    expect_is(p, "recordedplot")
    p[[1]][[6]][[2]]$col
})

test_that("plotInteractions(), which = 'network', metric = 'jaccard'", {
    p <- .(plotInteractions(x, "network", "jaccard"))
    expect_is(p, "recordedplot")
    p[[1]][[6]][[2]]$col
})

test_that("plotInteractions(), which = 'network', tunning", {
    p <- .(plotInteractions(
        x,
        which = "network",
        edge.color = "cyan",
        vertex.color = "pink",
        vertex.label.font = 2,
        vertex.label.color = "maroon"
    ))
    expect_is(p, "recordedplot")
    # Test edge color
    expect_equal(p[[1]][[6]][[2]]$col, "cyan")
    # Vertex label color and font
    expect_equal(p[[1]][[10]][[2]][[9]], "maroon")
    expect_equal(p[[1]][[10]][[2]][[10]], 2)
    # Vertex color
    expect_equal(p[[1]][[8]][[2]][[7]], "pink")
    p[[1]]
})
