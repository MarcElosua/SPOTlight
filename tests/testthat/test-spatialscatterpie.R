set.seed(321)
# Coordinates
x <- matrix(nrow = 10, data = c(1:10, 10:1))
rownames(x) <- paste0("spot", 1:nrow(x))
colnames(x) <- c("imagecol", "imagerow")
# Proportions
y <- replicate(m <- 5, runif(10, 0, 1))
y <- y / rowSums(y)
rownames(y) <- paste0("spot", 1:nrow(y))
colnames(y) <- paste0("type", 1:ncol(y))
# image
img <- paste0(system.file(package = "SPOTlight"), "/extdata/image.png")

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y)
    expect_equal(class(plt)[1], "gg")
})

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie - image", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y,
        img = img)
    expect_equal(class(plt)[1], "gg")
})

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie - type subset", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y,
        cell_types = colnames(y)[1:3])
    expect_equal(class(plt)[1], "gg")
})

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie - type subset", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y,
        scatterpie_alpha = 0.5)
    
    expect_equal(class(plt)[1], "gg")
    expect_lt(plt$layers[[1]]$aes_params$alpha, 1)
})

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie - type subset", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y,
        pie_scale = 0.1)
    expect_equal(class(plt)[1], "gg")
})

