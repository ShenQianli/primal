#' Plot function for S3 class "primal"
#'
#' Plot regularization path and parameter obtained from the algorithm.
#'
#' @param str An object with S3 class \code{"primal"}
#' @param n If \code{n = NULL}, three graph will be shown together. If \code{n} is a number, then the corresponding graph will be shown.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{Dantzig_solver}}, \code{\link{SparseSVM_solver}}
#' @export
plot.primal <- function(str, n = NULL, ...) {
    tt <- str$iterN
        opar <- par(no.readonly = TRUE)
    if (is.null(n)) {
        par(mfrow = c(1, 3), family = "serif")
        matplot(str$lambda[1:tt], t(str$beta[, 1:tt]), type = "l", 
                main = "Regularization Path", xlab = "Regularization Parameter", 
                ylab = "Coefficient", cex.main = 2, cex.lab = 1.6)
        matplot(1:tt, t(str$beta), type = "l", main = "Regularization Path", 
                xlab = "Iteration", ylab = "Coefficient", cex.main = 2, cex.lab = 1.6)
        plot(1:tt, str$lambda[1:tt], type = "l", main = "Value of Lambda along the Path", 
             xlab = "Iteration", ylab = "Lambda", cex.main = 2, cex.lab = 1.6)
        par(opar)
    } else {
        opar <- par(no.readonly = TRUE)
        par(family = "serif")
        switch(n, 
               matplot(str$lambda[1:tt], t(str$beta[, 1:tt]), type = "l", 
                       main = "Regularization Path", xlab = "Regularization Parameter", 
                       ylab = "Coefficient", cex.main = 2, cex.lab = 1.6), 
               matplot(1:tt, t(str$beta[, 1:tt]), type = "l", main = "Regularization Path", 
                       xlab = "Iteration", ylab = "Coefficient", cex.main = 2, cex.lab = 1.6),
               plot(1:tt, str$lambda[1:tt], type = "l", main = "Value of Lambda along the Path", 
                    xlab = "Iteration", ylab = "Lambda", cex.main = 2, cex.lab = 1.6)
               )
        par(opar)
        
    }
}


