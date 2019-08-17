#' Coef function for S3 class "primal"
#'
#' Print the estimated solution correspond to a specific parameter.
#'
#' @param str An object with S3 class \code{"primal"}.
#' @param n The index of the wanted parameter.
#' @seealso \code{\link{Dantzig_solver}}, \code{\link{SparseSVM_solver}}
#' @export
coef.primal <- function(str, n = NULL) {
    if (is.null(n))
        n <- str$iterN
    cat(str$type, " problem", "\n")
    cat("index: ", n)
    cat("\nlambda: ", str$lambda[n])
    cat("\ndegree of freedom: ", str$df[n])
    cat("\nbeta: \n")
    print(str$beta[, n])
    if(str$type == "SparseSVM"){
        cat("\nbeta0: ", str$beta0[n])
    }
}
