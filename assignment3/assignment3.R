#Load data
solar <- read.csv("assignment3/datasolar.csv")

# 2.1. 

#Parameters
mu <- 5.72 
ar1 <- 0.38
ar12 <- 0.94
ar13 <- - ar1 * ar12
e_var <- 0.22^2
k <- 12

#Introducing X_t
solar$xt <- log(solar$power) - mu


predict_power_with_ci <- function(X, ar1, ar12, ar13, mu, k=12, e_var) {
  X_pred <- numeric(k)

  #X_history keeps track of the last 13 values of X_t which are then updated in each iteration to include the new prediction
  X_history <- tail(X, 13)
  
  epsilon_t <- rnorm(k, mean = 0, sd = sqrt(e_var)) 
  
  # Standard deviation of the residuals (error variance)
  st_dev <- sqrt(e_var)

  # Critical value for 95% confidence interval (normal distribution)
  z_critical <- 1.96

  lower_bound <- numeric(k)
  upper_bound <- numeric(k)
  
  for (i in 1:k) {
    X_pred[i] <- ar1 * X_history[13] + ar12 * X_history[2] + ar13 * X_history[1] + epsilon_t[i]
    X_history <- c(X_history[-1], X_pred[i]) # it updates the history with the new prediction
    
    lower_bound[i] <- X_pred[i] - z_critical * st_dev
    upper_bound[i] <- X_pred[i] + z_critical * st_dev
  }
  
  #We convert predictions and bounds back to the original scale
  Y_pred <- exp(X_pred + mu)
  Y_lower <- exp(lower_bound + mu)
  Y_upper <- exp(upper_bound + mu)
  
  df <- data.frame(
    year = rep("2011", k),
    month = paste0(1:k),
    power = Y_pred,
    power_lower = Y_lower,
    power_upper = Y_upper,
    xt = X_pred
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


