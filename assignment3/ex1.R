library(astsa)
library(httpgd)
# Open server to view plots
hgd()

## 1 - Stability
# Functions
sim_ar2 <- function(phi1, phi2, seed = 2, n = 200, realizations = 5) {
  set.seed(seed)
  real_matrix <- matrix(NA, nrow = n, ncol = realizations)
  for (i in 1:realizations) {
    real_matrix[, i] <- arima.sim(
      n = n, model = list(order = c(2, 0, 0), ar = c(-phi1, -phi2))
    )
  }
  return(real_matrix)
}

plot_realiz <- function(real_matrix, phi1, phi2) {
  matplot(real_matrix, type = "l", lty = 1, col = seq_len(ncol(real_matrix)),
          main = paste("AR(2) with 5 realizations",
                       "for phi1 =", phi1, "and phi2 =", phi2),
          xlab = "Time", ylab = "Value")
  legend("topright", legend = paste("Realization", seq_len(ncol(real_matrix))),
         col = seq_len(ncol(real_matrix)), lty = 1, cex = 0.8)
}

plot_acf <- function(real_matrix, phi1, phi2, lags = 30) {
  emp_acfs <- matrix(NA, nrow = lags + 1, ncol = ncol(real_matrix))
  for (i in seq_len(ncol(real_matrix))) {
    emp_acfs[, i] <- acf(real_matrix[, i], lag.max = lags, plot = FALSE)$acf
  }

  emp_acf <- rowMeans(emp_acfs)
  theor_acf <- ARMAacf(ar = c(-phi1, -phi2), lag.max = lags)

  plot(emp_acf, type = "h", lty = 1, col = "blue",
       main = paste("Empirical vs Theoretical ACF",
                    "for phi1 =", phi1, "and phi2 =", phi2),
       xlab = "Lag", ylab = "ACF")

  offset <- 0:(length(theor_acf) - 1)
  segments(x0 = offset + 0.1, y0 = 0,
           x1 = offset + 0.1, y1 = theor_acf, col = "red", lwd = 2)

  legend("topright", legend = c("Empirical ACF", "Theoretical ACF"),
         col = c("blue", "red"), lty = c(1, 1), cex = 0.8)
}

# 1.1
phi1 <- -0.6
phi2 <- 0.5
real_matrix <- sim_ar2(phi1, phi2)
plot_realiz(real_matrix, phi1, phi2)

# 1.2
plot_acf(real_matrix, phi1, phi2)

# 1.3
phi1 <- -0.6
phi2 <- -0.3
plot_acf(real_matrix, phi1, phi2)

# 1.4
phi1 <- 0.6
phi2 <- -0.3
plot_acf(real_matrix, phi1, phi2)

# 1.5
phi1 <- -0.7
phi2 <- -0.3
plot_acf(real_matrix, phi1, phi2)

# 1.6
phi1 <- -0.75
phi2 <- -0.3
plot_acf(real_matrix, phi1, phi2)
