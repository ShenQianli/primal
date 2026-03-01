#' Plot a "primal" Object
#'
#' Plot the regularization path and parameters obtained from the parametric simplex method.
#'
#' @param x An object with S3 class \code{"primal"}.
#' @param n If \code{NULL}, all three plots are shown together. If \code{n} is a number (1, 2, or 3), only the corresponding plot is shown.
#' @param \dots Additional arguments (currently unused).
#' @seealso \code{\link{Dantzig_solver}}, \code{\link{SparseSVM_solver}}
#' @export
plot.primal <- function(x, n = NULL, ...) {
    tt <- x$iterN
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    if (is.null(n)) {
        par(mfrow = c(1, 3), family = "serif")
        matplot(x$lambda, t(x$beta), type = "l",
                main = "Regularization Path", xlab = "Regularization Parameter",
                ylab = "Coefficient", cex.main = 2, cex.lab = 1.6)
        matplot(1:tt, t(x$beta), type = "l", main = "Regularization Path",
                xlab = "Iteration", ylab = "Coefficient", cex.main = 2, cex.lab = 1.6)
        plot(1:tt, x$lambda, type = "l", main = "Value of Lambda along the Path",
             xlab = "Iteration", ylab = "Lambda", cex.main = 2, cex.lab = 1.6)
    } else {
        par(family = "serif")
        switch(n,
               matplot(x$lambda, t(x$beta), type = "l",
                       main = "Regularization Path", xlab = "Regularization Parameter",
                       ylab = "Coefficient", cex.main = 2, cex.lab = 1.6),
               matplot(1:tt, t(x$beta), type = "l", main = "Regularization Path",
                       xlab = "Iteration", ylab = "Coefficient", cex.main = 2, cex.lab = 1.6),
               plot(1:tt, x$lambda, type = "l", main = "Value of Lambda along the Path",
                    xlab = "Iteration", ylab = "Lambda", cex.main = 2, cex.lab = 1.6)
               )
    }
}
