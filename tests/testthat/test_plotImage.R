x <- "~/packages/SPOTlight/inst/extdata/SPOTlight.png"

test_that("plotImage,character", {
    # valid call
    y <- expect_silent(plotImage(x))
    expect_is(y, "ggplot")
    expect_equivalent(
        sapply(y$coordinates$limits, .subset, 2), 
        rev(dim(readPNG(x))[-3]))
    # invalid file type
    expect_error(plotImage(paste0(x, "x")))
    # invalid PNG
    x <- tempfile("image", fileext = ".png")
    cat("foo", file = x)
    expect_error(plotImage(x))
    file.remove(x)
})
