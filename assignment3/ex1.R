library(astsa)
library(ggplot2)
library(httpgd)
# Open server to view plots
hgd()

# 1. Stability

# Generate AR(2) realizations
sim_ar2 <- function(phi1, phi2, seed = 2, n = 200, realizations = 5) {
  set.seed(seed)
  # Create an empty matrix to store the realizations
  real_matrix <- matrix(NA, nrow = n, ncol = realizations)
  # Populate matrix with AR(2) realizations
  for (i in 1:realizations) {
    real_matrix[, i] <- arima.sim(
      n = n, model = list(order = c(2, 0, 0), ar = c(-phi1, -phi2)) # Sign is flipped accordig to implementation
    )
  }
  return(real_matrix)
}

# Plot realizations
plot_realiz <- function(real_matrix, phi1, phi2) {
  matplot(real_matrix, type = "l", lty = 1, col = seq_len(ncol(real_matrix)),
          main = paste("AR(2) with 5 realizations",
                       "for phi1 =", phi1, "and phi2 =", phi2),
          xlab = "Time", ylab = "Value")
  legend("topright", legend = paste("Realization", seq_len(ncol(real_matrix))),
         col = seq_len(ncol(real_matrix)), lty = 1, cex = 0.8)
}

# Plot empirical and theoretical ACF
plot_acf <- function(real_matrix, phi1, phi2, lags = 30) {
  # Get empirical ACFs of all 5 realizations
  emp_acfs <- matrix(NA, nrow = lags + 1, ncol = ncol(real_matrix))
  for (i in seq_len(ncol(real_matrix))) {
    emp_acfs[, i] <- acf(real_matrix[, i], lag.max = lags, plot = FALSE)$acf
  }

  # emp_acf <- rowMeans(emp_acfs)

  # Get theoretical ACF of the process
  theor_acf <- ARMAacf(ar = c(-phi1, -phi2), lag.max = lags)
  theor_acf_df <- as.data.frame(theor_acf)
  theor_acf_df$Lag <- 0:lags

  # Plot
  plot_df <- data.frame(Realization = character(), Lag = character(), ACF = character())
  for (col in 1:ncol(emp_acfs)) {
    rel_name <- paste('Realization', as.character(col), sep=" ")
    plot_df <- rbind(plot_df, data.frame(Realization = rel_name, Lag = 0:lags, ACF = emp_acfs[, col]))
  }

  ggplot() +
  geom_bar(data=plot_df, aes(fill=Realization, y=ACF, x=Lag), position="dodge", stat="identity", width=0.5) + 
  geom_line(data=theor_acf_df, aes(x=Lag, y=theor_acf, color="Theoretical ACF"), stat='Identity', color="black", size=1) +
  geom_point(data=theor_acf_df, aes(x=Lag, y=theor_acf), color="black", size=3) +
  guides(fill = guide_legend(override.aes = list(linetype = 0))) +
  ggtitle(paste("Theoretical vs Empirical ACFs",
                    "for phi1 =", phi1, "and phi2 =", phi2)) +
  theme_minimal() +
  theme(
        legend.title = element_blank(), 
        plot.title = element_text(face = "bold", size = (15), hjust = 0.5)
      )
  
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
