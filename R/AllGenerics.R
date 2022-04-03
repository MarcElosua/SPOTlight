#' @export
setGeneric(
    "SPOTlight", function(x, y, ...) 
        standardGeneric("SPOTlight"))

#' @export
setGeneric(
    "trainNMF", function(x, y, ...) 
        standardGeneric("trainNMF"))

#' @export
setGeneric(
    "runDeconvolution", function(x, ...) 
        standardGeneric("runDeconvolution"))

#' #' @export
#' setGeneric(
#'     "plotImage", function(x, ...) 
#'         standardGeneric("plotImage"))

#' @export
setGeneric(
    "plotTopicProfiles", function(x, y, ...) 
        standardGeneric("plotTopicProfiles"))

#' @export
setGeneric(
    "plotSpatialScatterpie", function(x, y, ...) 
        standardGeneric("plotSpatialScatterpie"))
