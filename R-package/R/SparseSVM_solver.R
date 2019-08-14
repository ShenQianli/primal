#-----------------------------------------------------------------------#
# Package: Parametric Simplex Method                                    #
# QuantileRegression_solver() : Solve given Sparse SVM problem in       #
# parametric simplex method                                             #
#-----------------------------------------------------------------------#

#' Solve given Sparse SVM problem in parametric simplex method
#'
#' See more details in \code{\link{psm}}
#' @param X \code{x} is an \code{n} by \code{d} data matrix
#' @param y \code{y} is a length \code{n} response vector
#' @param max_it This is the number of the maximum path length one would like to achieve. The default length is \code{50}.
#' @param lambda_threshold The parametric simplex method will stop when the calculated parameter is smaller than lambda. The default value is \code{0.01}.
#' @export
SparseSVM_solver <- function(X, y, max_it = 50, lambda_threshold = 0.01){

  n <- nrow(X)
  d <- ncol(X)
  t <- 0
  lambdalist <- rep(0, max_it)
  x_list <- matrix(0, d, max_it)
  y_list <- rep(0, max_it)
  x0 <- rep(0, max_it)

  str <- .C("SparseSVM_api", as.integer(n), as.integer(d),
            as.double(X), as.double(y),
            as.integer(max_it), as.double(lambda_threshold),
            as.integer(t), as.double(lambdalist),
            as.double(x_list), as.double(y_list),
            as.double(x0), PACKAGE = "PRIMAL")
  t <- unlist(str[7])
  x_list <- matrix(unlist(str[9], d, t))

  return(list(
    type = "SparseSVM",
    BETA = x_list,
    BETA0 = unlist(str[11]),
    value = unlist(str[10]),
    iterN_given = max_it,
    iterN = t,
    lambdalist = unlist(str[8])
    )
  )

}
