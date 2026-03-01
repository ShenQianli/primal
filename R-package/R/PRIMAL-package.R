#' @useDynLib PRIMAL, .registration = TRUE
#' @importFrom stats runif
#' @importFrom graphics par image plot lines matplot
#' @importFrom Matrix Matrix t
NULL

#' Parametric Simplex Method for Sparse Learning
#'
#' Implements the parametric simplex method for a variety of sparse
#' learning problems including Dantzig selector, sparse SVM, compressed
#' sensing, and quantile regression.
#'
#' \tabular{ll}{
#'   Package: \tab PRIMAL\cr
#'   Type: \tab Package\cr
#'   Version: \tab 1.1.0\cr
#'   Date: \tab 2026-02-27\cr
#' }
#' The package provides four main solvers:\cr
#' (1) Dantzig selector: \code{\link{Dantzig_solver}}\cr
#' (2) Sparse SVM: \code{\link{SparseSVM_solver}}\cr
#' (3) Compressed sensing: \code{\link{CompressedSensing_solver}}\cr
#' (4) Quantile regression: \code{\link{QuantileRegression_solver}}\cr
#' @docType package
#' @aliases primal-package
#' @author Qianli Shen, Zichong Li, Tuo Zhao \cr
#' @seealso \code{\link{plot.primal}}, \code{\link{print.primal}}, \code{\link{coef.primal}}
"_PACKAGE"
#> [1] "_PACKAGE"