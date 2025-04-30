source("assignment4/functions/kalmanfilter.R")
source("assignment4/functions/myLogLikFun.R")

#1.1
# Set parameters
a <- 0.9
b <- 1
sigma <- 1
X0 <- 5
n <- 100
num_realizations <- 5

# Set seed for reproducibility
set.seed(123)

# Create a matrix to store the results
# Each column will represent one realization
simulations <- matrix(0, nrow = n+1, ncol = num_realizations)
simulations
# Set initial values (X0) for all realizations
simulations[1,] <- X0

# Simulate the AR(1) process for each realization
for (j in 1:num_realizations) {
  for (t in 1:n) {
    # Generate random error from normal distribution
    error <- rnorm(1, mean = 0, sd = sigma)
    
    # Calculate the next value in the sequence using the AR(1) formula
    simulations[t+1, j] <- a * simulations[t, j] + b + error
  }
}

# Create a data frame for easier plotting
sim_df <- data.frame(time = 0:n)
for (j in 1:num_realizations) {
  sim_df[paste0("realization_", j)] <- simulations[, j]
}

# Plot the realizations
png("AR1_realizations.png", width = 10, height = 6, units = "in", res = 300)

# Your existing plot code
plot(0:n, simulations[, 1], type = "l", col = 1, 
     ylim = range(simulations), 
     xlab = "Time", ylab = "Value", 
     main = "5 Realizations of AR(1) Process")
for (j in 2:num_realizations) {
  lines(0:n, simulations[, j], col = j)
}
legend("bottomright", legend = paste0("Series ", 1:num_realizations), 
       col = 1:num_realizations, lty = 1)

# Close the device to save the file
dev.off()

#1.2
# Set parameters
a <- 0.9
b <- 1
sigma1 <- 1  # Process noise standard deviation
sigma2 <- 1  # Observation noise standard deviation
X0 <- 5      # Initial value
n <- 100     # Number of time steps

# Create vectors to store the results
X <- numeric(n+1)  # True state (Xt)
Y <- numeric(n+1)  # Observations (Yt)

# Set initial value
X[1] <- X0
Y[1] <- X0 + rnorm(1, mean = 0, sd = sigma2)  # First observation with noise

# Simulate the AR(1) process and observations
for (t in 1:n) {
  # Generate process noise
  e1 <- rnorm(1, mean = 0, sd = sigma1)
  
  # Calculate the next state value
  X[t+1] <- a * X[t] + b + e1
  
  # Generate observation noise
  e2 <- rnorm(1, mean = 0, sd = sigma2)
  
  # Calculate the observation
  Y[t+1] <- X[t+1] + e2
}

# Create time vector
time <- 0:n

# Plot the true state and the observations
png("state_and_observations.png", width = 10, height = 6, units = "in", res = 300)

plot(time, X, type = "l", col = "blue", lwd = 2,
     ylim = range(c(X, Y)),  # Set y-axis range to include both X and Y
     xlab = "Time", ylab = "Value",
     main = "True State (X) and Noisy Observations (Y)")

# Add observations
lines(time, Y, type = "p", col = "red", pch = 16, cex = 0.7)

# Add legend
legend("topright", 
       legend = c("True State (X)", "Observations (Y)"),
       col = c("blue", "red"), 
       lty = c(1, NA), 
       pch = c(NA, 16),
       lwd = c(2, NA),
       cex = 0.8)

dev.off()

# Also display in the current graphics device
plot(time, X, type = "l", col = "blue", lwd = 2,
     ylim = range(c(X, Y)),
     xlab = "Time", ylab = "Value",
     main = "True State (X) and Noisy Observations (Y)")
lines(time, Y, type = "p", col = "red", pch = 16, cex = 0.7)
legend("topright", 
       legend = c("True State (X)", "Observations (Y)"),
       col = c("blue", "red"), 
       lty = c(1, NA), 
       pch = c(NA, 16),
       lwd = c(2, NA),
       cex = 0.8)

#1.3
# Now apply the Kalman filter
theta <- c(a, b, sigma1)  # [0.9, 1, 1]
R <- sigma2^2            # Observation noise variance is sigma2^2 = 1

# Apply the Kalman filter
kf_results <- myKalmanFilter(
  y = Y[-1],           # Remove the first observation (at t=0)
  theta = theta,       # Parameters [a, b, sigma1]
  R = R,               # Measurement noise variance
  x_prior = X0,        # Initial state
  P_prior = 1        # Initial variance (small value as we know X0 exactly)
)

# Create the plot
time <- 0:n
time_kf <- 1:n  # Kalman filter times (starts at t=1)

# Calculate the 95% confidence intervals (±1.96 standard deviations)
upper_ci <- kf_results$x_pred + 1.96 * sqrt(kf_results$P_pred)
lower_ci <- kf_results$x_pred - 1.96 * sqrt(kf_results$P_pred)

# Plot
png("kalman_filter_results.png", width = 10, height = 6, units = "in", res = 300)
plot(time, X, type = "l", col = "blue", lwd = 2,
     ylim = range(c(X, Y, upper_ci, lower_ci)),
     xlab = "Time", ylab = "Value",
     main = "Kalman Filter Results")

# Add observations
points(time, Y, col = "red", pch = 16, cex = 0.7)

# Add predicted state
lines(time_kf, kf_results$x_pred, col = "green", lwd = 2)

# Add confidence intervals
lines(time_kf, upper_ci, col = "green", lty = 2)
lines(time_kf, lower_ci, col = "green", lty = 2)

# Add a semi-transparent green area for the confidence interval
polygon(c(time_kf, rev(time_kf)), c(lower_ci, rev(upper_ci)),
        col = rgb(0, 1, 0, 0.2), border = NA)

# Add legend
legend("topright", 
       legend = c("True State (X)", "Observations (Y)", "Predicted State", "95% CI"),
       col = c("blue", "red", "green", "green"), 
       lty = c(1, NA, 1, 2), 
       pch = c(NA, 16, NA, NA),
       lwd = c(2, NA, 2, 1),
       cex = 0.8)
dev.off()

#1.4

set.seed(123)  # For reproducibility

estimate_ar1_parameters <- function(a, b, sigma1, sigma2 = 1, X0 = 5, n = 100, num_realizations = 100, 
                                    plot_title_prefix = "", save_plots = FALSE, plot_filename_prefix = "ar1_params") {
  
  # Create a list to store parameter estimates
  param_estimates <- list(
    a = numeric(num_realizations),
    b = numeric(num_realizations),
    sigma1 = numeric(num_realizations)
  )
  
  # Function to simulate a single realization of the AR(1) process with observations
  simulate_ar1_with_obs <- function(a, b, sigma1, sigma2, X0, n) {
    X <- numeric(n+1)  # True state (Xt)
    Y <- numeric(n+1)  # Observations (Yt)
    
    # Set initial value
    X[1] <- X0
    Y[1] <- X0 + rnorm(1, mean = 0, sd = sigma2)
    
    # Simulate the AR(1) process and observations
    for (t in 1:n) {
      # Generate process noise
      e1 <- rnorm(1, mean = 0, sd = sigma1)
      
      # Calculate the next state value
      X[t+1] <- a * X[t] + b + e1
      
      # Generate observation noise
      e2 <- rnorm(1, mean = 0, sd = sigma2)
      
      # Calculate the observation
      Y[t+1] <- X[t+1] + e2
    }
    
    return(list(X = X, Y = Y))
  }
  
  # Loop over each realization
  for (i in 1:num_realizations) {
    # Simulate a realization
    sim_data <- simulate_ar1_with_obs(a, b, sigma1, sigma2, X0, n)
    Y <- sim_data$Y
    
    # Initial parameter guess
    initial_theta <- c(0.5, 0.5, 0.5)  # Starting with neutral values
    
    # Optimize the negative log-likelihood
    opt_result <- optim(
      par = initial_theta,
      fn = myLogLikFun,
      y = Y[-1],  # Remove the first observation (at t=0)
      R = sigma2^2,  # Observation noise variance
      x_prior = X0,  # Initial state
      P_prior = 1,   # Initial variance (small as we know X0 exactly)
      method = "L-BFGS-B",  # Optimization method
      lower = c(0, -10, 0.01),  # Lower bounds for parameters
      upper = c(10, 10, 5)    # Upper bounds for parameters
    )
    
    # Store parameter estimates
    param_estimates$a[i] <- opt_result$par[1]
    param_estimates$b[i] <- opt_result$par[2]
    param_estimates$sigma1[i] <- opt_result$par[3]
    
    # Print progress
    if (i %% 10 == 0) {
      cat("Completed realization", i, "of", num_realizations, "\n")
    }
  }
  
  # Convert parameter estimates to a data frame for easier plotting
  param_df <- data.frame(
    a = param_estimates$a,
    b = param_estimates$b,
    sigma1 = param_estimates$sigma1
  )
  
  # Create a summary of the parameter estimates
  summary_stats <- data.frame(
    Parameter = c("a", "b", "sigma1"),
    True_Value = c(a, b, sigma1),
    Mean = c(mean(param_df$a), mean(param_df$b), mean(param_df$sigma1)),
    Median = c(median(param_df$a), median(param_df$b), median(param_df$sigma1)),
    SD = c(sd(param_df$a), sd(param_df$b), sd(param_df$sigma1)),
    Min = c(min(param_df$a), min(param_df$b), min(param_df$sigma1)),
    Max = c(max(param_df$a), max(param_df$b), max(param_df$sigma1))
  )
  
  # Create and save boxplots if requested
  if (save_plots) {
    png(paste0(plot_filename_prefix, "_individual.png"), width = 12, height = 4, units = "in", res = 300)
  }
  
  # Create boxplots
  par(mfrow = c(1, 3))
  
  # Boxplot for a
  boxplot(param_df$a, main = paste0(plot_title_prefix, " Estimates of a"), ylab = "Value")
  abline(h = a, col = "red", lwd = 2)
  text(1.3, a + 0.05, "True value", col = "red")
  
  # Boxplot for b
  boxplot(param_df$b, main = paste0(plot_title_prefix, " Estimates of b"), ylab = "Value")
  abline(h = b, col = "red", lwd = 2)
  text(1.3, b + 0.1, "True value", col = "red")
  
  # Boxplot for sigma1
  boxplot(param_df$sigma1, main = paste0(plot_title_prefix, " Estimates of sigma1"), ylab = "Value")
  abline(h = sigma1, col = "red", lwd = 2)
  text(1.3, sigma1 + 0.1, "True value", col = "red")
  
  if (save_plots) {
    dev.off()
    png(paste0(plot_filename_prefix, "_combined.png"), width = 8, height = 6, units = "in", res = 300)
  }
  
  if (save_plots) {
    dev.off()
  }
  
  # Return the results
  return(list(
    summary_stats = summary_stats,
    param_df = param_df
  ))
}

results1 <- estimate_ar1_parameters(
  a = 1, 
  b = 0.9, 
  sigma1 = 1,
  plot_title_prefix = "Case 1:",
  save_plots = TRUE,
  plot_filename_prefix = "case1_params"
)

results2 <- estimate_ar1_parameters(
  a = 5, 
  b = 0.9, 
  sigma1 = 1,
  plot_title_prefix = "Case 2:",
  save_plots = TRUE,
  plot_filename_prefix = "case2_params"
)

# Try a third parameter set
results3 <- estimate_ar1_parameters(
  a = 1, 
  b = 0.9, 
  sigma1 = 5,
  plot_title_prefix = "Case 3:",
  save_plots = TRUE,
  plot_filename_prefix = "case3_params"
)

#1.5
set.seed(123)  # For reproducibility
# Function to simulate AR(1) process with t-distributed noise
simulate_ar1_with_t_noise <- function(a, b, sigma1, sigma2, X0, n, df) {
  X <- numeric(n+1)  # True state (Xt)
  Y <- numeric(n+1)  # Observations (Yt)
  
  # Set initial value
  X[1] <- X0
  Y[1] <- X0 + rnorm(1, mean = 0, sd = sigma2)
  
  # Simulate the AR(1) process and observations
  for (t in 1:n) {
    # Generate process noise from t-distribution
    e1 <- sigma1 * rt(1, df = df)
    
    # Calculate the next state value
    X[t+1] <- a * X[t] + b + e1
    
    # Generate observation noise (still Gaussian)
    e2 <- rnorm(1, mean = 0, sd = sigma2)
    
    # Calculate the observation
    Y[t+1] <- X[t+1] + e2
  }
  
  return(list(X = X, Y = Y))
}

# Set parameters
a <- 1
b <- 0.9
sigma1 <- 1
sigma2 <- 1
X0 <- 5
n <- 100
num_simulations <- 100
df_values <- c(100, 5, 2, 1)  # Degrees of freedom to test

# Create plot to compare t and normal distributions
png("t_vs_normal.png", width = 10, height = 6, units = "in", res = 300)

# Define x-range for plotting
x <- seq(-5, 5, length.out = 1000)

# Plot standard normal distribution first
plot(x, dnorm(x, mean = 0, sd = 1), type = "l", col = "black", lwd = 2,
     xlab = "x", ylab = "Density", 
     main = "Comparison of Normal and t Distributions",
     ylim = c(0, 0.45))  # Setting appropriate y-axis limits

# Define degrees of freedom values
df_values <- c(100, 5, 2, 1)
colors <- c("blue", "red", "green", "purple")

# Add t-distributions
for (i in 1:length(df_values)) {
  lines(x, dt(x, df = df_values[i]), col = colors[i], lwd = 2)
}

# Add legend
legend("topright", 
       legend = c("Normal", paste("t(", df_values, ")", sep = "")),
       col = c("black", colors), 
       lwd = 2)

dev.off()

# Store simulations results
all_simulations <- list()

# Run simulations for each df value
for (i in 1:length(df_values)) {
  df <- df_values[i]
  simulations <- list()
  
  for (j in 1:num_simulations) {
    sim_data <- simulate_ar1_with_t_noise(a, b, sigma1, sigma2, X0, n, df)
    simulations[[j]] <- sim_data
  }
  
  all_simulations[[i]] <- simulations
}


estimate_ar1_parameters <- function(a, b, sigma1, sigma2 = 1, X0 = 5, n = 100,df, num_realizations = 100, 
                                    plot_title_prefix = "", save_plots = FALSE, plot_filename_prefix = "ar1_params") {
  
  # Create a list to store parameter estimates
  param_estimates <- list(
    a = numeric(num_realizations),
    b = numeric(num_realizations),
    sigma1 = numeric(num_realizations)
  )
  
  # Loop over each realization
  for (i in 1:num_realizations) {
    # Simulate a realization
    sim_data <- simulate_ar1_with_t_noise(a, b, sigma1, sigma2, X0, n, df)
    Y <- sim_data$Y
    
    # Initial parameter guess
    initial_theta <- c(0.5, 0.5, 0.5)  # Starting with neutral values
    
    # Optimize the negative log-likelihood
    opt_result <- optim(
      par = initial_theta,
      fn = myLogLikFun,
      y = Y[-1],  # Remove the first observation (at t=0)
      R = sigma2^2,  # Observation noise variance
      x_prior = X0,  # Initial state
      P_prior = 1,   # Initial variance (small as we know X0 exactly)
      method = "L-BFGS-B",  # Optimization method
      lower = c(0, -10, 0.01),  # Lower bounds for parameters
      upper = c(10, 10, 5)    # Upper bounds for parameters
    )
    
    # Store parameter estimates
    param_estimates$a[i] <- opt_result$par[1]
    param_estimates$b[i] <- opt_result$par[2]
    param_estimates$sigma1[i] <- opt_result$par[3]
    
    # Print progress
    if (i %% 10 == 0) {
      cat("Completed realization", i, "of", num_realizations, "\n")
    }
  }
  
  # Convert parameter estimates to a data frame for easier plotting
  param_df <- data.frame(
    a = param_estimates$a,
    b = param_estimates$b,
    sigma1 = param_estimates$sigma1
  )
  
  # Create a summary of the parameter estimates
  summary_stats <- data.frame(
    Parameter = c("a", "b", "sigma1"),
    True_Value = c(a, b, sigma1),
    Mean = c(mean(param_df$a), mean(param_df$b), mean(param_df$sigma1)),
    Median = c(median(param_df$a), median(param_df$b), median(param_df$sigma1)),
    SD = c(sd(param_df$a), sd(param_df$b), sd(param_df$sigma1)),
    Min = c(min(param_df$a), min(param_df$b), min(param_df$sigma1)),
    Max = c(max(param_df$a), max(param_df$b), max(param_df$sigma1))
  )
  
  # Create and save boxplots if requested
  if (save_plots) {
    png(paste0(plot_filename_prefix, "_individual.png"), width = 12, height = 4, units = "in", res = 300)
  }
  
  # Create boxplots
  par(mfrow = c(1, 3))
  
  # Boxplot for a
  boxplot(param_df$a, main = paste0(plot_title_prefix, " Estimates of a"), ylab = "Value")
  abline(h = a, col = "red", lwd = 2)
  text(1.3, a + 0.05, "True value", col = "red")
  
  # Boxplot for b
  boxplot(param_df$b, main = paste0(plot_title_prefix, " Estimates of b"), ylab = "Value")
  abline(h = b, col = "red", lwd = 2)
  text(1.3, b + 0.1, "True value", col = "red")
  
  # Boxplot for sigma1
  boxplot(param_df$sigma1, main = paste0(plot_title_prefix, " Estimates of sigma1"), ylab = "Value")
  abline(h = sigma1, col = "red", lwd = 2)
  text(1.3, sigma1 + 0.1, "True value", col = "red")
  
  if (save_plots) {
    dev.off()
    png(paste0(plot_filename_prefix, "_combined.png"), width = 8, height = 6, units = "in", res = 300)
  }
  
  if (save_plots) {
    dev.off()
  }
  
  # Return the results
  return(list(
    summary_stats = summary_stats,
    param_df = param_df
  ))
}
for (i in 1:length(df_values)) {
  df <- df_values[i]
  results <- estimate_ar1_parameters(
    a = 1, 
    b = 0.9, 
    sigma1 = 1,
    df = df,
    plot_title_prefix = paste("Case with t-distributed noise (df =", df, "):"),
    save_plots = TRUE,
    plot_filename_prefix = paste("case_t_params", df, sep = "_")
  )}
