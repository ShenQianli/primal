if (!require("flare")) install.packages("flare")
if (!require("PRIMAL")) install.packages("PRIMAL")
if (!require("fastclime")) install.packages("fastclime")

library(flare)
library(PRIMAL)
library(fastclime)


#################### test dantzig selector time######################
cat("**************Test dantzig selector time****************")
fastclime_time <- rep(0, 6)
PRIMAL_time <- rep(0, 6)
for(i in 1:7){
    data <- dantzig.generator(n = 200, d = 1000*i, sparsity = 0.03)
    fastclime_time[i] <- system.time(dantzig(data$X0, data$y, nlambda = i*10000, lambda = 2 *
    sqrt(log(i*1000)/200)))
    ans <- Dantzig_solver(data$X0, data$y, max_it = i*10000, lambda_threshold = 2 *
    sqrt(log(i*1000)/200))
    PRIMAL_time[i] <-ans$runtime
}

time_data <- matrix(c(seq(1000, 7000, 1000), fastclime_time, PRIMAL_time), ncol = 3)
opar <- par(no.readonly = TRUE)
par(family = "serif")
plot(time_data[, 1:2], type = "l", col = "orange", lwd = 3, main = "Dantzig selector", xlab = "Dimension", ylab = "CPU Time(s)", cex.main = 2, cex.lab = 1.5)
legend("topleft", inset = 0.05, c("PRIMAL", "fastclime"), col = c("blue", "orange"), lwd = c(3, 3), cex = 1.5)
lines(time_data[, c(1, 3)], type = "l", lwd = 3, col = "blue")
par(opar)


#################### test compressed sensing time######################
cat("**************Test compressed sensing time****************")

# function to solve compressed sensing problem using fastclime package #
test_fastclime_c <- function(X, y, lambda) {
    c_obj <- rep(-1, 2 * ncol(X))
    A <- rbind(cbind(X, -X), cbind(-X, X))
    b_rhs <- rbind(data$y, -data$y)
    c_bar <- rep(0, 2 * ncol(X))
    b_bar <- rep(1, 2 * nrow(X))
    opt <- paralp(obj = c_obj, mat = A, rhs = b_rhs, obj_bar = c_bar, 
                  rhs_bar = b_bar, lambda = lambda)
    return(opt)
}

fastclime_time2 <- rep(0, 9)
PRIMAL_time2 <- rep(0, 9)
flare_time2 <- rep(0, 9)

for (i in 1:9) {
    for (j in 1:trialN) {
        data <- dantzig.generator(n = 200, d = 100 * (i + 1), sparsity = 0.1)
        fastclime_time2[i] <- fastclime_time2[i] + system.time(test_fastclime_c(data$X0, data$y, lambda = 2 * sqrt(log(100 *
            (i + 1))/200)))
        
        PRIMAL_time2[i] <- PRIMAL_time2[i] + system.time(ans <- CompressedSensing_solver(data$X0, data$y, max_it = 10000, 
            lambda_threshold = 2 * sqrt(log(100 * (i + 1))/200)))
        
    }
}

# plot out the time result#
time_data2 <- matrix(c(seq(200, 1000, 100), fastclime_time2/trialN, PRIMAL_time2/trialN), ncol = 3)
opar <- par(no.readonly = TRUE)
par(family = "serif")
plot(time_data2[, 1:2], type = "l", col = "orange", lwd = 3, main = "Compressed Sensing", xlab = "Dimension", ylab = "CPU Time(s)", 
    ylim = c(0, max(fastclime_time2/trialN) + 0.1), cex.main = 2, cex.lab = 1.5)
legend("topleft", inset = 0.05, c("PRIMAL", "fastclime"), col = c("blue", "orange"), lwd = c(3, 3), cex = 1.5)
lines(time_data2[, c(1, 3)], type = "l", lwd = 3, col = "blue")
par(opar)




# plot out the time result together#
opar <- par(no.readonly = TRUE)
par(family = "serif", mfrow = c(1, 2))
plot(time_data[, 1:2], type = "l", col = "orange", lwd = 3, main = "Dantzig selector", xlab = "Dimension", ylab = "CPU Time(s)", cex.main = 2, cex.lab = 1.5)
legend("topleft", inset = 0.05, c("PRIMAL", "fastclime"), col = c("blue", "orange"), lwd = c(3, 3,
    3), cex = 1.5)
lines(time_data[, c(1, 3)], type = "l", lwd = 3, col = "blue")
lines(time_data[, c(1, 4)], type = "l", lwd = 3, col = "brown")
plot(time_data2[, 1:2], type = "l", col = "orange", lwd = 3, main = "Compressed Sensing", xlab = "Dimension", ylab = "CPU Time(s)", 
    ylim = c(0, max(fastclime_time2/trialN) + 0.1), cex.main = 2, cex.lab = 1.5)
legend("topleft", inset = 0.05, c("PRIMAL", "fastclime"), col = c("blue", "orange"), lwd = c(3, 3), cex = 1.5)
lines(time_data2[, c(1, 3)], type = "l", lwd = 3, col = "blue")
par(opar)

