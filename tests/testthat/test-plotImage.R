x_path <- paste0(system.file(package = "SPOTlight"), "/extdata/image.png")
# plotImage() ----
test_that("plotImage path", {
    # image
    x <- x_path
    plt <- plotImage(x = x)
    
    expect_equal(class(plt)[1], "gg")
})


# plotImage() ----
test_that("plotImage array", {
    # image
    x <- png::readPNG(x_path)
    plt <- plotImage(x = x)
    
    expect_equal(class(plt)[1], "gg")
})

test_that("plotImage Seurat", {
    # image
    if (! "stxBrain.SeuratData" %in% SeuratData::InstalledData()$Dataset) {
        SeuratData::InstallData("stxBrain.SeuratData")
    }
    
    x <- SeuratData::LoadData("stxBrain.SeuratData")
    
    plt <- plotImage(x = x)
    
    expect_equal(class(plt)[1], "gg")
})

test_that("plotImage SPE", {
    # image
    library(ExperimentHub)
    eh <- ExperimentHub()        # initialize hub instance
    q <- query(eh, "TENxVisium") # retrieve 'TENxVisiumData' records
    id <- q$ah_id[1]             # specify dataset ID to load
    x <- eh[[id]]
        
    plt <- plotImage(x = x)
    
    expect_equal(class(plt)[1], "gg")
})
