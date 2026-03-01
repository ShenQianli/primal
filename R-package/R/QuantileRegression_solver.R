#-----------------------------------------------------------------------#
# Package: PaRametric sImplex Method for spArse Learning                #
# QuantileRegression_solver() : Solve given quantile regression problem #
# in parametric simplex method                                          #
#-----------------------------------------------------------------------#

#' Solve a Quantile Regression Problem Using the Parametric Simplex Method
#'
#' @param X An \code{n} by \code{d} data matrix.
#' @param y A length \code{n} response vector.
#' @param max_it Maximum number of iterations for the solution path. The default is \code{50}.
#' @param lambda_threshold The algorithm stops when the regularization parameter falls below this threshold. The default value is \code{0.01}.
#' @param tau The quantile level for the regression. The default is \code{0.5}.
#' @note The returned value is similar to that of \code{\link{Dantzig_solver}}.
#' @seealso \code{\link{primal-package}}, \code{\link{Dantzig_solver}}
#' @export
QuantileRegression_solver <- function(X, y, max_it = 50, lambda_threshold = 0.01, tau = 0.5) {
    begt <- Sys.time()
    n <- nrow(X)
    d <- ncol(X)
    iter_count <- 0
    lambdalist <- rep(0, max_it)
    x_list <- matrix(0, d, max_it)
    y_list <- rep(0, max_it)
    res <- .C("QuantileRegression_api", n = as.integer(n), d = as.integer(d),
              X = as.double(t(X)), y = as.double(y), tau = as.double(tau),
              max_it = as.integer(max_it), lambda_threshold = as.double(lambda_threshold),
              iter_count = as.integer(iter_count),
              lambdalist = as.double(lambdalist),
              x_list = as.double(x_list), y_list = as.double(y_list), PACKAGE = "PRIMAL")
    iter_count <- res$iter_count
    x_list <- matrix(res$x_list[1:(d * iter_count)], d, iter_count)
    df <- colSums(x_list != 0)
    runt <- Sys.time() - begt
    ans <- list(type = "Quantile Regression",
                data = X,
                response = y,
                beta = Matrix(x_list),
                df = df,
                value = res$y_list[1:iter_count],
                iterN = iter_count,
                lambda = res$lambdalist[1:iter_count],
                runtime = runt)
    class(ans) <- "primal"
    return(ans)
}
