#-----------------------------------------------------------------------#
# Package: PaRametric sImplex Method for spArse Learning                # 
# QuantileRegression_solver() : Solve given quantile regression problem # 
# in parametric simplex method                                          #
#-----------------------------------------------------------------------#

#' Solve given quantile regression problem in parametric simplex method
#'
#' @param X \code{x} is an \code{n} by \code{d} data matrix
#' @param y \code{y} is a length \code{n} response vector
#' @param max_it This is the number of the maximum path length one would like to achieve. The default length is \code{50}.
#' @param lambda_threshold The parametric simplex method will stop when the calculated parameter is smaller than lambda. The default value is \code{0.01}.
#' @param tau The quantile number you want. The default quantile is \code{0.5}
#' @note The returned value is similar with \code{Dantzig_solver()}
#' @seealso \code{\link{primal-package}}, \code{\link{Dantzig_solver}}
#' @export
QuantileRegression_solver <- function(X, y, max_it = 50, lambda_threshold = 0.01, tau = 0.5) {
    begt <- Sys.time()
    n <- nrow(X)
    d <- ncol(X)
    t <- 0
    lambdalist <- rep(0, max_it)
    x_list <- matrix(0, d, max_it)
    y_list <- rep(0, max_it)
    str <- .C("QuantileRegression_api", as.integer(n), as.integer(d), as.double(t(X)), 
              as.double(y), as.double(tau), as.integer(max_it), as.double(lambda_threshold), 
              as.integer(t), as.double(lambdalist), as.double(x_list), as.double(y_list), PACKAGE = "PRIMAL")
    t <- unlist(str[8])
    x_list <- matrix(unlist(str[10])[1:(d * t)], d, t)
    df <- c()
    for (i in 1:t) {
        df[i] = sum(x_list[, i] != 0)
    }
    runt <- Sys.time() - begt
    ans <- list(type = "Quantile Regression", 
                data = X,
                response = y,
                beta = x_list, 
                df = df, 
                value = unlist(str[11])[1:t], 
                iterN = t, 
                lambda = unlist(str[9])[1:t], 
                runtime = runt)
    class(ans) <- "primal"
    return(ans)
}
