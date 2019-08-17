#' @useDynLib PRIMAL, .registration = TRUE
#' @importFrom stats runif
#' @importFrom graphics par image plot lines matplot
NULL

#' Parametric Simplex Method for Sparse Learning
#'
#' A package for parametric simplex method for sparse learning
#'
#' \tabular{ll}{
#'   Package: \tab PRIMAL\cr
#'   Type: \tab Package\cr
#'   Version: \tab 1.0.0\cr
#'   Date: \tab 2019-08-15\cr
#' }
#' The package "PRIMAL" provides 4 main functions:\cr
#' (1) The dantzig selector solver applying simplex method. Please refer to \code{\link{Dantzig_solver}}.\cr
#' (2) The sparse SVM solver applying simplex method. Please refer to \code{\link{SparseSVM_solver}}.\cr
#' (3) The compressed sensing solver. Please refer to \code{\link{CompressedSensing_solver}}.\cr
#' (4) The quantile regression solver. Please refer to \code{\link{QuantileRegression_solver}}.\cr
#' @docType package
#' @aliases primal-package
#' @author Qianli Shen, Zichong Li \cr
#' @seealso \code{\link{plot.primal}}, \code{\link{print.primal}}, \code{\link{coef.primal}}
"_PACKAGE"
#> [1] "_PACKAGE"