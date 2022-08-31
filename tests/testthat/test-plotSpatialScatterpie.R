set.seed(321)
# Coordinates
x <- matrix(nrow = 10, data = c(seq_len(10), 10:1))
rownames(x) <- paste0("spot", seq_len(nrow(x)))
colnames(x) <- c("coord_x", "coord_y")
# Proportions
y <- replicate(m <- 5, runif(10, 0, 1))
y <- y / rowSums(y)
rownames(y) <- paste0("spot", seq_len(nrow(y)))
colnames(y) <- paste0("type", seq_len(ncol(y)))
# image
img <- paste0(system.file(package = "SPOTlight"), "/extdata/SPOTlight.png")

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie with matrix and bad colnames", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y
    )
    expect_equal(class(plt)[1], "gg")
})

colnames(x) <- c("imagecol", "imagerow")
# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie with matrix and bad colnames", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y
    )
    expect_equal(class(plt)[1], "gg")
})

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie - image", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y,
        img = img
    )
    expect_equal(class(plt)[1], "gg")
})

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie - type subset", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y,
        cell_types = colnames(y)[seq_len(3)]
    )
    expect_equal(class(plt)[1], "gg")
})

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie - alpha", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y,
        scatterpie_alpha = 0.5
    )

    expect_equal(class(plt)[1], "gg")
    expect_lt(plt$layers[[1]]$aes_params$alpha, 1)
})

# plotSpatialScatterpie() ----
test_that("plotSpatialScatterpie - pie_scale", {
    plt <- plotSpatialScatterpie(
        x = x,
        y = y,
        pie_scale = 0.1
    )
    expect_equal(class(plt)[1], "gg")
})

library(SpatialExperiment)
example(read10xVisium, echo = FALSE)
# Proportions
spe_y <- replicate(m <- 5, runif(ncol(spe), 0, 1))
spe_y <- spe_y / rowSums(spe_y)
rownames(spe_y) <- colnames(spe)
colnames(spe_y) <- paste0("type", seq_len(ncol(spe_y)))

# plotSpatialScatterpie() img TRUE----
test_that("plotSpatialScatterpie - image", {
    plt <- plotSpatialScatterpie(
        x = spe,
        y = spe_y,
        img = TRUE
    )
    expect_equal(class(plt)[1], "gg")
    # Make sure there is an image
    expect_true(is(plt$layers[[1]]$geom, "GeomCustomAnn"))
})

# plotSpatialScatterpie() img TRUE----
test_that("plotSpatialScatterpie - spots on image", {
    plt <- plotSpatialScatterpie(
        x = spe,
        y = spe_y,
        img = TRUE
    )
    expect_equal(class(plt)[1], "gg")
    # Make sure there is an image
    expect_true(is(plt$layers[[1]]$geom, "GeomCustomAnn"))
    
    # Check the spots are on within the image coordinates
    x_y_min_max <- plt$layers[[1]]$geom_params
    point_df <- plt$layers[[2]]$data
    expect_true(max(point_df$coord_x) <= x_y_min_max$xmax)
    expect_true(min(point_df$coord_x) >= x_y_min_max$xmin)
    expect_true(max(point_df$coord_y) <= x_y_min_max$ymax)
    expect_true(min(point_df$coord_y) >= x_y_min_max$ymin)
})


# plotSpatialScatterpie() img TRUE----
test_that("plotSpatialScatterpie - rotate 90 degrees", {
    plt <- plotSpatialScatterpie(
      x = spe,
      y = spe_y,
      img = TRUE,
      degrees = 90
    )
    expect_equal(class(plt)[1], "gg")
    # Make sure there is an image
    expect_true(is(plt$layers[[1]]$geom, "GeomCustomAnn"))
    
    # Check the spots are on within the image coordinates
    x_y_min_max <- plt$layers[[1]]$geom_params
    point_df <- plt$layers[[2]]$data
    expect_true(max(point_df$coord_x) <= x_y_min_max$xmax)
    expect_true(min(point_df$coord_x) >= x_y_min_max$xmin)
    expect_true(max(point_df$coord_y) <= x_y_min_max$ymax)
    expect_true(min(point_df$coord_y) >= x_y_min_max$ymin)
})

# plotSpatialScatterpie() img TRUE----
test_that("plotSpatialScatterpie - mirror veritcal", {
  plt <- plotSpatialScatterpie(
    x = spe,
    y = spe_y,
    img = TRUE,
    axis = "v"
  )
  expect_equal(class(plt)[1], "gg")
  # Make sure there is an image
  expect_true(is(plt$layers[[1]]$geom, "GeomCustomAnn"))
  
  # Check the spots are on within the image coordinates
  x_y_min_max <- plt$layers[[1]]$geom_params
  point_df <- plt$layers[[2]]$data
  expect_true(max(point_df$coord_x) <= x_y_min_max$xmax)
  expect_true(min(point_df$coord_x) >= x_y_min_max$xmin)
  expect_true(max(point_df$coord_y) <= x_y_min_max$ymax)
  expect_true(min(point_df$coord_y) >= x_y_min_max$ymin)
})

test_that("plotSpatialScatterpie - mirror horizontal", {
  plt <- plotSpatialScatterpie(
    x = spe,
    y = spe_y,
    img = TRUE,
    axis = "h"
  )
  expect_equal(class(plt)[1], "gg")
  # Make sure there is an image
  expect_true(is(plt$layers[[1]]$geom, "GeomCustomAnn"))
  
  # Check the spots are on within the image coordinates
  x_y_min_max <- plt$layers[[1]]$geom_params
  point_df <- plt$layers[[2]]$data
  expect_true(max(point_df$coord_x) <= x_y_min_max$xmax)
  expect_true(min(point_df$coord_x) >= x_y_min_max$xmin)
  expect_true(max(point_df$coord_y) <= x_y_min_max$ymax)
  expect_true(min(point_df$coord_y) >= x_y_min_max$ymin)
})
