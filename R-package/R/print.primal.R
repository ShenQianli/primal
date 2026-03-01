#' Print a "primal" Object
#'
#' Print summary information including the model type, regularization parameter path, and degrees of freedom along the solution path.
#'
#' @param x An object with S3 class \code{"primal"}.
#' @param \dots Additional arguments (currently unused).
#' @seealso \code{\link{Dantzig_solver}}, \code{\link{SparseSVM_solver}}
#' @export
print.primal <- function(x, ...) {
    cat("\n*******Parametric Simplex Method solving ")
    cat(x$type, "problem*********\n")
    cat("iteration times = ", x$iterN, "\n")
    cat("lambda list:\n")
    print(signif(x$lambda, digits = 5))
    cat("Degree of freedom:", x$df[1], "----->", x$df[x$iterN], "\n")
    unit <- units.difftime(x$runtime)
    cat("Runtime:", x$runtime, " ", unit, "\n")
}
