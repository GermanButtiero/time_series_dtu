#Load data
solar <- read.csv("assignment3/datasolar.csv")

mu <- 5.72 
ar1 <- -0.38  
ar12 <- -0.94  
ar13 <- ar1 * ar12  
e_var <- 0.22^2
k <- 12

# 2.1
solar$xt <- log(solar$power) - mu
# Calculate one-step ahead prediction errors (residuals)
# We need to use the formula:
# ε̂(t+1|t) = X(t+1) + φ1*X(t) + Φ1*X(t-11) + φ1*Φ1*X(t-12)

residuals <- NA

for(t in 14:nrow(solar)) {
      residuals[t] <- solar$xt[t] + 
                       ar1 * solar$xt[t-1] + 
                       ar12 * solar$xt[t-12] + 
                       ar13 * solar$xt[t-13]
}

# Model validation - check i.i.d. assumption
# 1. Basic statistics
mean_residuals <- mean(residuals, na.rm = TRUE)
var_residuals <- var(residuals, na.rm = TRUE)
cat("Mean of residuals:", mean_residuals, "\n")
cat("Variance of residuals:", var_residuals, "\n")
cat("Expected variance:", e_var, "\n")

# 2. Check independence (autocorrelation)
png("acf_residuals.png", width=800, height=600)
acf(residuals, na.action = na.pass, main = "ACF of Residuals")
dev.off()

png("histogram_residuals.png", width=800, height=600)
hist(residuals, breaks=30, main="Histogram of Residuals", xlab="Residuals", col="lightblue", border="black")
dev.off()

png("scatter_residuals.png", width=800, height=600)
plot(residuals, type="p", pch=16, col="blue", main="Residuals Scatter Plot", xlab="Time", ylab="Residuals")
abline(h=0, col="red", lty=2)
dev.off()

png("qqplot_residuals.png", width=800, height=600)
qqnorm(residuals, main="Q-Q Plot of Residuals")
qqline(residuals, col="red", lwd=2)
dev.off()

# Shapiro-Wilk test
shapiro_test <- shapiro.test(residuals)
cat("Shapiro-Wilk test p-value:", shapiro_test$p.value, "\n")
#since the p-value is 0.1515613 > 0.05, this means that we cannot reject the null hypothesis
#that the residuals are normally distributed.

predict_power_with_ci <- function(X, ar1, ar12, ar13, mu, k, e_var) {
  X_pred <- numeric(k)
  X_history <- tail(X, 13)  

  # Calculate prediction variance for each horizon
  var_k <- numeric(k)
  # For AR(1) model, we need to calculate the coefficients for the prediction variance
  # psi_1 = ar1, psi_2 = ar1^2, ..., psi_(k-1) = ar1^(k-1)
  psi_coef <- numeric(k-1)
  if (k > 1) {
    psi_coef[1] <- ar1  
    for (j in 2:(k-1)) {
      psi_coef[j] <- ar1 * psi_coef[j-1]
    }
  }
  
  for (i in 1:k) {
    X_pred[i] <- -ar1 * X_history[13] - ar12 * X_history[2] - ar13 * X_history[1]

    X_history <- c(X_history[-1], X_pred[i])
    
    # variance using equation 5.151
    # For AR(1): var_k[i] = e_var * (1 + psi_1^2 + psi_2^2 + ... + psi_(i-1)^2)
    if (i == 1) {
      var_k[i] <- e_var
    } else {
      var_k[i] <- e_var * (1 + sum(psi_coef[1:(i-1)]^2))
    }
  }
  
  # Critical value for 95% confidence interval
  z_critical <- qnorm(0.975)  
  
  # Calculate confidence intervals
  lower_bound <- X_pred - z_critical * sqrt(var_k)
  upper_bound <- X_pred + z_critical * sqrt(var_k)
  
  # Transform back to power (need to account for log transformation properly)
  Y_pred <- exp(X_pred + mu)
  
  # For log-transformed data, the correct CI transformation is:
  Y_lower <- exp(X_pred + mu - z_critical * sqrt(var_k))
  Y_upper <- exp(X_pred + mu + z_critical * sqrt(var_k))
  
  # Create dataframe for results
  df <- data.frame(
    year = rep(2011, k),
    month = 1:k,
    power = Y_pred,
    xt = X_pred,
    power_lower = Y_lower,
    power_upper = Y_upper
  )
  
  return(df)
}

set.seed(1)
predictions <- predict_power_with_ci(solar$xt, ar1, ar12, ar13, mu, k, e_var)

solar$power_lower <- NA #to match predictions structure
solar$power_upper <- NA #to match predictions structure

df <- rbind(solar, predictions)

plotit <- function(df, filename="plot2.1.png") {
  png(filename, width=800, height=600)
  
  n_obs <- nrow(df) - k  
  
  layout(rbind(1, 2:3))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7,0))
  
  
  plot(df$power, type="l", col="blue", pch=16, ylab="Power", xlab="Time", main="Observed vs Forecasted Power")
  lines((n_obs+1):nrow(df), df$power[(n_obs+1):nrow(df)], col="red", pch=16)
  polygon(c(1:nrow(df), rev(1:nrow(df))),
          c(df$power_lower, rev(df$power_upper)),
          col=rgb(0, 0, 1, alpha=0.2), border=NA) 
  
  legend("topleft", legend=c("Observed", "Forecasted", "95% CI"), col=c("blue", "red", rgb(0, 0, 1, alpha=0.2)), pch=c(16, 16, NA), lty=c(NA, NA, 1), lwd=c(2, 2, 2))

  # ACF and PACF plots
  acf(df$power, lag.max=50, lwd=2, main="ACF of Power")
  pacf(df$power, lag.max=50, lwd=2, main="PACF of Power")
  dev.off()
}
plotit(df, "plot3.2.png")

#2.4
#Comment: would you trust the forecast? Do you think the prediction intervals have correct
#Since the AR(1) model only accounts for the immediate past (time t-1) and does not consider 
#long-term patterns like seasonality, we cannot fully trust the forecast, especially for 
#longer horizons. The plot supports this, as the prediction intervals widen significantly 
#after the first few observations, reflecting the increasing uncertainty in the forecast as
# we move further from the observed data.