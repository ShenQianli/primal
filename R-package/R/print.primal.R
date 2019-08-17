#' Print function for S3 class "primal"
#'
#' Print the information about the model usage, the parameter path, degree of freedom of the solution path.
#'
#' @param str An object with S3 class \code{"primal"}.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{Dantzig_solver}}, \code{\link{SparseSVM_solver}}
#' @export
print.primal <- function(str, ...) {
    cat("\n*******Parametric Simplex Method solving ")
    cat(str$type, "problem*********\n")
    cat("iteration times = ", str$iterN, "\n")
    cat("lambda list:\n")
    print(signif(str$lambda, digits = 5))
    cat("Degree of freedom:", str$df[1], "----->", str$df[str$iterN], "\n")
    if (units.difftime(str$runtime) == "secs") 
        unit <- "secs"
    if (units.difftime(str$runtime) == "mins") 
        unit <- "mins"
    if (units.difftime(str$runtime) == "hours") 
        unit <- "hours"
    cat("Runtime:", str$runtime, " ", unit, "\n")
}

