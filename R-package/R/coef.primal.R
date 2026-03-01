#' Extract Coefficients from a "primal" Object
#'
#' Extract and display the estimated coefficients corresponding to a specific regularization parameter.
#'
#' @param object An object with S3 class \code{"primal"}.
#' @param n The index of the regularization parameter along the solution path. If \code{NULL}, the last iteration is used.
#' @param \dots Additional arguments (currently unused).
#' @seealso \code{\link{Dantzig_solver}}, \code{\link{SparseSVM_solver}}
#' @export
coef.primal <- function(object, n = NULL, ...) {
    if (is.null(n))
        n <- object$iterN
    cat(object$type, " problem", "\n")
    cat("index: ", n)
    cat("\nlambda: ", object$lambda[n])
    cat("\ndegree of freedom: ", object$df[n])
    cat("\nbeta: \n")
    print(object$beta[, n])
    if(object$type == "SparseSVM"){
        cat("\nbeta0: ", object$beta0[n])
    }
}
