#' @importFrom matrixStats colMedians
#' @importFrom NMF coef
topicProfiles <- function(mod, groups) {
    df <- data.frame(t(coef(mod)))
    dfs <- split(df, groups)
    res <- vapply(
        dfs, function(df)
            colMedians(as.matrix(df)),
        numeric(ncol(df))
    )
    rownames(res) <- names(dfs)
    return(t(res))
}
