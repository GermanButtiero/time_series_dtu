source("assignment4/functions/kalmanfilter.R")
source("assignment4/functions/myLogLikFun.R")

#1.1 
set.seed(123)
a <- 0.9
b <- 1
sigma <- 1
X0 <- 5
n <- 100
num_realizations <- 5

simulations <- matrix(0, nrow = n+1, ncol = num_realizations)
simulations[1,] <- X0

for (j in 1:num_realizations) {
  for (t in 1:n) {
    error <- rnorm(1, mean = 0, sd = sigma)
    simulations[t+1, j] <- a * simulations[t, j] + b + error
  }
}

sim_df <- data.frame(time = 0:n)
for (j in 1:num_realizations) {
  sim_df[paste0("realization_", j)] <- simulations[, j]
}

png("assignment4/plots/AR1_realizations.png", width = 10, height = 6, units = "in", res = 300)

plot(0:n, simulations[, 1], type = "l", col = 1, 
     ylim = range(simulations), 
     xlab = "Time", ylab = "Value", 
     main = "5 Realizations of AR(1) Process")
for (j in 2:num_realizations) {
  lines(0:n, simulations[, j], col = j)
}
legend("bottomright", legend = paste0("Series ", 1:num_realizations), 
       col = 1:num_realizations, lty = 1)
dev.off()

#1.2
set.seed(123) 
a <- 0.9
b <- 1
sigma1 <- 1  
sigma2 <- 1 
X0 <- 5      
n <- 100  

X <- numeric(n+1)  
Y <- numeric(n+1)  

X[1] <- X0
Y[1] <- X0 + rnorm(1, mean = 0, sd = sigma2)  

#
for (t in 1:n) {
  e1 <- rnorm(1, mean = 0, sd = sigma1)
  X[t+1] <- a * X[t] + b + e1

  e2 <- rnorm(1, mean = 0, sd = sigma2)
  Y[t+1] <- X[t+1] + e2
}

time <- 0:n

png("assignment4/plots/state_and_observations.png", width = 10, height = 6, units = "in", res = 300)

plot(time, X, type = "l", col = "blue", lwd = 2,
     ylim = range(c(X, Y)),  # Set y-axis range to include both X and Y
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
dev.off()

#1.3
set.seed(123)  
theta <- c(a, b, sigma1) 
R <- sigma2^2            


kf_results <- myKalmanFilter(
  y = Y[-1],           
  theta = theta,       
  R = R,               
  x_prior = X0,       
  P_prior = 1   
)

time <- 0:n
time_kf <- 1:n  

upper_ci <- kf_results$x_pred + 1.96 * sqrt(kf_results$P_pred)
lower_ci <- kf_results$x_pred - 1.96 * sqrt(kf_results$P_pred)

png("assignment4/plots/kalman_filter_results.png", width = 10, height = 6, units = "in", res = 300)
plot(time, X, type = "l", col = "blue", lwd = 2,
     ylim = range(c(X, Y, upper_ci, lower_ci)),
     xlab = "Time", ylab = "Value",
     main = "Kalman Filter Results")
points(time, Y, col = "red", pch = 16, cex = 0.7)
lines(time_kf, kf_results$x_pred, col = "green", lwd = 2)
lines(time_kf, upper_ci, col = "green", lty = 2)
lines(time_kf, lower_ci, col = "green", lty = 2)
polygon(c(time_kf, rev(time_kf)), c(lower_ci, rev(upper_ci)),
        col = rgb(0, 1, 0, 0.2), border = NA)
legend("topright", 
       legend = c("True State (X)", "Observations (Y)", "Predicted State", "95% CI"),
       col = c("blue", "red", "green", "green"), 
       lty = c(1, NA, 1, 2), 
       pch = c(NA, 16, NA, NA),
       lwd = c(2, NA, 2, 1),
       cex = 0.8)
dev.off()

#1.4
set.seed(123)  

estimate_ar1_parameters <- function(a, b, sigma1, sigma2 = 1, X0 = 5, n = 100, num_realizations = 100, 
                                    plot_title_prefix = "", save_plots = FALSE, plot_filename_prefix = "ar1_params") {
  
  param_estimates <- list(
    a = numeric(num_realizations),
    b = numeric(num_realizations),
    sigma1 = numeric(num_realizations)
  )
  simulate_ar1_with_obs <- function(a, b, sigma1, sigma2, X0, n) {
    X <- numeric(n+1)  
    Y <- numeric(n+1)  
    
    X[1] <- X0
    Y[1] <- X0 + rnorm(1, mean = 0, sd = sigma2)
    
    for (t in 1:n) {
      e1 <- rnorm(1, mean = 0, sd = sigma1)
      X[t+1] <- a * X[t] + b + e1
      
      e2 <- rnorm(1, mean = 0, sd = sigma2)
      Y[t+1] <- X[t+1] + e2
    }
    
    return(list(X = X, Y = Y))
  }
  
  for (i in 1:num_realizations) {
    sim_data <- simulate_ar1_with_obs(a, b, sigma1, sigma2, X0, n)
    Y <- sim_data$Y
    
    initial_theta <- c(0.5, 0.5, 0.5)  
    
    opt_result <- optim(
      par = initial_theta,
      fn = myLogLikFun,
      y = Y[-1],  
      R = sigma2^2,  
      x_prior = X0, 
      P_prior = 1,  
      method = "L-BFGS-B", 
      lower = c(0, -10, 0.01),  
      upper = c(10, 10, 5) 
    )
    
    param_estimates$a[i] <- opt_result$par[1]
    param_estimates$b[i] <- opt_result$par[2]
    param_estimates$sigma1[i] <- opt_result$par[3]
    
    }
  
  param_df <- data.frame(
    a = param_estimates$a,
    b = param_estimates$b,
    sigma1 = param_estimates$sigma1
  )
  
  summary_stats <- data.frame(
    Parameter = c("a", "b", "sigma1"),
    True_Value = c(a, b, sigma1),
    Mean = c(mean(param_df$a), mean(param_df$b), mean(param_df$sigma1)),
    Median = c(median(param_df$a), median(param_df$b), median(param_df$sigma1)),
    SD = c(sd(param_df$a), sd(param_df$b), sd(param_df$sigma1)),
    Min = c(min(param_df$a), min(param_df$b), min(param_df$sigma1)),
    Max = c(max(param_df$a), max(param_df$b), max(param_df$sigma1))
  )
  
  if (save_plots) {
    png(paste0("assignment4/plots/",plot_filename_prefix, "_individual.png"), width = 12, height = 4, units = "in", res = 300)
  }
  
  par(mfrow = c(1, 3))
  boxplot(param_df$a, main = paste0(plot_title_prefix, " Estimates of a"), ylab = "Value")
  abline(h = a, col = "red", lwd = 2)
  text(1.3, a + 0.05, "True value", col = "red")
  
  boxplot(param_df$b, main = paste0(plot_title_prefix, " Estimates of b"), ylab = "Value")
  abline(h = b, col = "red", lwd = 2)
  text(1.3, b + 0.1, "True value", col = "red")
  
  boxplot(param_df$sigma1, main = paste0(plot_title_prefix, " Estimates of sigma1"), ylab = "Value")
  abline(h = sigma1, col = "red", lwd = 2)
  text(1.3, sigma1 + 0.1, "True value", col = "red")
  
  if (save_plots) {
    dev.off()
  }
  
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

results3 <- estimate_ar1_parameters(
  a = 1, 
  b = 0.9, 
  sigma1 = 5,
  plot_title_prefix = "Case 3:",
  save_plots = TRUE,
  plot_filename_prefix = "case3_params"
)

#1.5
set.seed(123)  

simulate_ar1_with_t_noise <- function(a, b, sigma1, sigma2, X0, n, df) {
  X <- numeric(n+1)  
  Y <- numeric(n+1)  
  
  X[1] <- X0
  Y[1] <- X0 + rnorm(1, mean = 0, sd = sigma2)
  
  for (t in 1:n) {
    e1 <- sigma1 * rt(1, df = df)
    X[t+1] <- a * X[t] + b + e1
    
    e2 <- rnorm(1, mean = 0, sd = sigma2)
    Y[t+1] <- X[t+1] + e2
  }
  
  return(list(X = X, Y = Y))
}

a <- 1
b <- 0.9
sigma1 <- 1
sigma2 <- 1
X0 <- 5
n <- 100
num_simulations <- 100
df_values <- c(100, 5, 2, 1)  

png("assignment4/plots/t_vs_normal.png", width = 10, height = 6, units = "in", res = 300)

x <- seq(-5, 5, length.out = 1000)

plot(x, dnorm(x, mean = 0, sd = 1), type = "l", col = "black", lwd = 2,
     xlab = "x", ylab = "Density", 
     main = "Comparison of Normal and t Distributions",
     ylim = c(0, 0.45))  

df_values <- c(100, 5, 2, 1)
colors <- c("blue", "red", "green", "purple")

for (i in 1:length(df_values)) {
  lines(x, dt(x, df = df_values[i]), col = colors[i], lwd = 2)
}
legend("topright", 
       legend = c("Normal", paste("t(", df_values, ")", sep = "")),
       col = c("black", colors), 
       lwd = 2)
dev.off()

all_simulations <- list()

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
  
  param_estimates <- list(
    a = numeric(num_realizations),
    b = numeric(num_realizations),
    sigma1 = numeric(num_realizations)
  )
  
  for (i in 1:num_realizations) {
    sim_data <- simulate_ar1_with_t_noise(a, b, sigma1, sigma2, X0, n, df)
    Y <- sim_data$Y
    
    initial_theta <- c(0.5, 0.5, 0.5) 
    
    opt_result <- optim(
      par = initial_theta,
      fn = myLogLikFun,
      y = Y[-1],  
      R = sigma2^2,  
      x_prior = X0, 
      P_prior = 1, 
      method = "L-BFGS-B", 
      lower = c(0, -10, 0.01), 
      upper = c(10, 10, 5) 
    )
    
    param_estimates$a[i] <- opt_result$par[1]
    param_estimates$b[i] <- opt_result$par[2]
    param_estimates$sigma1[i] <- opt_result$par[3]

  }
  
  param_df <- data.frame(
    a = param_estimates$a,
    b = param_estimates$b,
    sigma1 = param_estimates$sigma1
  )
  
  summary_stats <- data.frame(
    Parameter = c("a", "b", "sigma1"),
    True_Value = c(a, b, sigma1),
    Mean = c(mean(param_df$a), mean(param_df$b), mean(param_df$sigma1)),
    Median = c(median(param_df$a), median(param_df$b), median(param_df$sigma1)),
    SD = c(sd(param_df$a), sd(param_df$b), sd(param_df$sigma1)),
    Min = c(min(param_df$a), min(param_df$b), min(param_df$sigma1)),
    Max = c(max(param_df$a), max(param_df$b), max(param_df$sigma1))
  )
  
  if (save_plots) {
    png(paste0("assignment4/plots/", plot_filename_prefix, "_individual.png"), width = 12, height = 4, units = "in", res = 300)
  }
  
  par(mfrow = c(1, 3))

  boxplot(param_df$a, main = paste0(plot_title_prefix, " Estimates of a"), ylab = "Value")
  abline(h = a, col = "red", lwd = 2)
  text(1.3, a + 0.05, "True value", col = "red")
  
  boxplot(param_df$b, main = paste0(plot_title_prefix, " Estimates of b"), ylab = "Value")
  abline(h = b, col = "red", lwd = 2)
  text(1.3, b + 0.1, "True value", col = "red")
  
  boxplot(param_df$sigma1, main = paste0(plot_title_prefix, " Estimates of sigma1"), ylab = "Value")
  abline(h = sigma1, col = "red", lwd = 2)
  text(1.3, sigma1 + 0.1, "True value", col = "red")
  
  if (save_plots) {
    dev.off()
  }
  
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
