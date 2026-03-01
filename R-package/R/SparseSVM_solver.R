#-----------------------------------------------------------------------#
# Package: PaRametric sImplex Method for spArse Learning                #
# SparseSVM_solver() : Solve given Sparse SVM problem in                #
# parametric simplex method                                             #
#-----------------------------------------------------------------------#

#' Solve a Sparse SVM Problem Using the Parametric Simplex Method
#'
#' @param X An \code{n} by \code{d} data matrix.
#' @param y A length \code{n} response vector with values in \{-1, 1\}.
#' @param max_it Maximum number of iterations for the solution path. The default is \code{50}.
#' @param lambda_threshold The algorithm stops when the regularization parameter falls below this threshold. The default value is \code{0.01}.
#' @return
#' An object with S3 class \code{"primal"} is returned:
#' \item{data}{
#'   The \code{n} by \code{d} data matrix from the input.
#' }
#' \item{response}{
#'   The length \code{n} response vector from the input.
#' }
#' \item{beta}{
#'   A matrix of coefficient estimates whose columns correspond to regularization parameters along the solution path.
#' }
#' \item{beta0}{
#'   A vector of intercept estimates corresponding to regularization parameters along the solution path.
#' }
#' \item{df}{
#'   The degrees of freedom (number of nonzero coefficients) along the solution path.
#' }
#' \item{value}{
#'   The sequence of optimal objective function values corresponding to each \code{lambda}.
#' }
#' \item{iterN}{
#'   The number of iterations performed.
#' }
#' \item{lambda}{
#'   The sequence of regularization parameters obtained along the solution path.
#' }
#' \item{type}{
#'   The type of the problem, e.g., \code{"Dantzig"} or \code{"SparseSVM"}.
#' }
#' @seealso \code{\link{primal-package}}
#' @export
SparseSVM_solver <- function(X, y, max_it = 50, lambda_threshold = 0.01) {
    begt <- Sys.time()
    n <- nrow(X)
    d <- ncol(X)
    iter_count <- 0
    lambdalist <- rep(0, max_it)
    x_list <- matrix(0, d, max_it)
    y_list <- rep(0, max_it)
    x0 <- rep(0, max_it)

    res <- .C("SparseSVM_api", n = as.integer(n), d = as.integer(d),
              X = as.double(t(X)), y = as.double(y),
              max_it = as.integer(max_it), lambda_threshold = as.double(lambda_threshold),
              iter_count = as.integer(iter_count),
              lambdalist = as.double(lambdalist),
              x_list = as.double(x_list), y_list = as.double(y_list),
              x0 = as.double(x0), PACKAGE = "PRIMAL")
    iter_count <- res$iter_count
    x_list <- matrix(res$x_list[1:(d * iter_count)], d, iter_count)
    df <- colSums(x_list != 0)
    runt <- Sys.time() - begt
    ans <- list(type = "SparseSVM",
                data = X,
                response = y,
                beta = Matrix(x_list),
                beta0 = res$x0[1:iter_count],
                df = df,
                value = res$y_list[1:iter_count],
                iterN = iter_count,
                lambda = res$lambdalist[1:iter_count],
                runtime = runt)
    class(ans) <- "primal"
    return(ans)
}
