#install.packages("dplyr")  # Install (only needed once)
library(dplyr)  # Load the package
library(tidyverse)
library(ggplot2)
library(httpgd)
hgd()

#getting data
Dtrain
Dtest

#2.2
#plot of the train data with the estimated mean
plot(Dtrain$year, Dtrain$total)

fitlm <- lm(total ~ year, data = Dtrain)
summary(fitlm)
abline(fitlm, col = "red")

# values of the estimated parameters and standard errors
theta_hat <- coef(fitlm)
theta_hat
print(paste("The estimated value of the parameters are: ", theta_hat))

theta_hat_se <- summary(fitlm)$coef[, 2]
theta_hat_se
print(paste("The estimated value of the st errors are: ", theta_hat_se))

#2.3
# Make a forecast with intervals
Dtest$forecast <- predict(fitlm, newdata = Dtest, interval = "prediction", level=0.95)
Dtest$forecast #create table from this

#2.4 plot everything with prediction intervals
ggplot(Dtrain, aes(x = time, y = total)) +
  geom_point(col = "red") +
  geom_line(aes(y = predict(fitlm)), col = "red", size = 0.8) +
  geom_point(data = Dtest, aes(x = time, y = forecast[,1]), col = "blue", size = 0.5) +
  geom_ribbon(data = Dtest, aes(x = time, ymin = forecast[,2], ymax = forecast[,3]), inherit.aes=FALSE,fill = "red", alpha = 0.2) +
  geom_point(data=Dtest, aes(x = time, y = total), col = "blue", size = 0.8) # add the total for the test data - used in 2.5

#2.5: Is it a good forecast? Comments:
# - the forecast seems to overshoot the actual values, as the forecast is higher than the actual values
# - the high slope does not match the actual development of the data over time, especially in the period of the test data
# - most of the forecasted values are outside of the prediction interval, which means that the forecast is not very accurate

# supplementary analysis - display the actual and forecasted values in a table
comparison_table <- data.frame(
  Actual = Dtest$total,
  Forecasted = Dtest$forecast[,1]
)
print(comparison_table)

#2.6
e_ols <- Dtest$total - Dtest$forecast[,1]
e_ols

qqnorm(e_ols)
qqline(e_ols)
# the residuals fall close to the reference line in the middle of the plot. 
# However, they are significantly far from the reference line in the beginning and the end of the plot.

#plot residuals versus year
ggplot(Dtest, aes(x = time, y = e_ols)) +
  geom_point(col = "blue") +
  geom_line(aes(y = e_ols), col = "blue", size = 0.5)

##it doesnt seem like it matches the assumptions, as the errors should be distributed around the mean (0)
# In contrast to the assumption, the errors on the plot are negative due to the forecast overshooting the actual values.


# Part 3 - Weighted least squares

## Define variables
n <- length(Dtrain$time) # number of observations
p <- 2 # number of parameters
X <- cbind(1, Dtrain$year) # design matrix
y <- Dtrain$total # response vector
X_test <- cbind(1, Dtest$year) # design matrix for test data
y_test <- Dtest$total # response vector for test data
lambda <- 0.9 # exponential weights


# 3.1 -- check slides
## OLS -> identity matrix of size 72 x 72
## WLS -> diagonal matrix of size 72 x 72 with 1/(lambda^(n-1)) on the diagonal)

# 3.2
## Calculate weights
weights <- lambda^((n-1):0)
## Plot the weights
barplot(weights, names=1:72)

# 3.3
print(sum(weights)) # for WLS
print(72) # for OLS -> sum of weights is equal to the number of observations

# 3.4
SIGMA <- diag(n)
diag(SIGMA) <- 1/weights
## Estimate parameters with WLS
theta_WLS <- solve(t(X)%*%solve(SIGMA)%*%X)%*%(t(X)%*%solve(SIGMA)%*%y)
print(theta_WLS)

# 3.5
## Train estimates from WLS
yhat_wls <- X%*%theta_WLS

## Test predictions from WLS
y_pred_wls <- X_test%*%theta_WLS

## Prediction intervals
e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%solve(SIGMA)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - p))
Vmatrix_pred <- sigma2_wls * (1 + (X_test %*% solve(t(X)%*%solve(SIGMA)%*%X)) %*% t(X_test) )
y_pred_lwr_wls <- y_pred_wls - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
y_pred_upr_wls <- y_pred_wls + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))

## Add variables to dataframe
Dtrain$yhat_wls <- yhat_wls
Dtest$y_pred_wls <- y_pred_wls
Dtest$y_pred_lwr_wls <- y_pred_lwr_wls
Dtest$y_pred_upr_wls <- y_pred_upr_wls

### Get estimates from OLS
SIGMA_ols <- diag(n)

## Estimate parameters with WLS
theta_OLS <- solve(t(X)%*%solve(SIGMA_ols)%*%X)%*%(t(X)%*%solve(SIGMA_ols)%*%y)
yhat_ols <- X%*%theta_OLS

## Test predictions from WLS
y_pred_ols <- X_test%*%theta_OLS

## Add to dataframe
Dtrain$yhat_ols <- yhat_ols


## Plot data and predictions
# plot WLS (blue) together with OLS (red) and true test data (black):
ggplot(Dtrain, aes(x=year, y=total)) +
  geom_point(col="red") + 

  geom_line(aes(y=yhat_ols), col="red", size=.5, linetype=2) +
  geom_point(data=Dtest, aes(x=year,y=forecast[,1]), col="red", size=.5) +
  geom_ribbon(data=Dtest, aes(x=year,ymin=forecast[,2], ymax=forecast[,3]), inherit.aes=FALSE, alpha=0.1, fill="red") +

  geom_point(data=Dtest, aes(x=year,y=total), col="black") +

  geom_line(aes(y=yhat_wls), col="blue", size=.5) +
  geom_point(data=Dtest, aes(x=year,y=y_pred_wls), col="blue", size=.5) +
  geom_ribbon(data=Dtest, aes(x=year,ymin=y_pred_lwr_wls, ymax=y_pred_upr_wls), inherit.aes=FALSE, alpha=0.2, fill="blue") +
  
  coord_cartesian(xlim = c(2018,2025))


# Part 4 - RLS and optimization of lambda

# 4.2
D <- Dtrain
X <- cbind(1, D$year) 
y <- D$total 
n <-length(X[,1]) # number of observations

# Initialize R0
# We do this in order to make the matrix invertible
p <- 2  
R_0 <- diag(0.1, p, p)

# Function to implement Recursive Least Squares (RLS)
RLS_estimate <- function(X, y, R_0) {
  # Number of observations and parameters
  n <- length(y)
  p <- ncol(X)  

  # Initialize Theta matrix to store estimates
  Theta <- matrix(NA, nrow=n, ncol=p)

  # Initialize R and h
  R <- R_0   #it allows us to make the matrix invertible
  h <- rep(0, p)  # Initial h is a vector of zeros because theta_0 = 0

  # Loop
  for (i in 1:n) {
    x <- X[i,]
    y_i <- y[i]

    # Update R and h
    R <- R + x %*% t(x)  # R_t = R_{t-1} + x_t * x_t^T
    h <- h + x * y_i  # h_t = h_{t-1} + x_t * y_t

    # Estimate parameters 
    Theta[i,] <- solve(R) %*% h  # θ_t = R_t^-1 * h_t
  }

  return(Theta)
}
Theta_RLS <- RLS_estimate(X, y, R_0)
Theta_RLS[1:3,]

# Do you think it is intuitive to understand the details in the
#matrix calculations? If yes, give a short explanaition


# 4.3

#difference between estimates of OLS and RLS
print(paste("The difference between the estimates of OLS and RLS is: ", Theta_RLS[n,1] - theta_hat[1]))

#The values are very far apart from each other. This is caused by the initial value of R_0, which is set to a very low value.
#This causes the estimates to be very different from the OLS estimates. A way to reduce this difference is to reduce the initial value of R_0. 
#A larger R_0 will cause the estimates to be more influenced by the initial value, while a smaller R_0 will cause the estimates to be more influenced by the data.

R_0 <- diag(0.001, p, p)
Theta_RLS <- RLS_estimate(X, y, R_0)
print(paste("The difference between the estimates of OLS and RLS is: ", Theta_RLS[n,1] - theta_hat[1]))

#4.4 
# Function to implement Recursive Least Squares (RLS) with lambda
RLS_estimate_with_forgetting <- function(X, y, R_0, lambda) {
  n <- length(y)  
  p <- ncol(X)    
  
  # Initialize Theta, Rt, and h_t
  Theta <- matrix(NA, nrow=n, ncol=p)
  R <- R_0
  h_t <- rep(0, p)  
  
  # Loop through each time step
  for (t in 1:n) {
    x_t <- X[t, ]
    y_t <- y[t]
    
    # Update h_t recursively
    h_t <- h_t + x_t * y_t
    
    # Update R_t with forgetting factor
    R <- lambda * R + x_t %*% t(x_t)
    
    # Calculate theta_t (parameter estimates)
    Theta[t, ] <- solve(R) %*% h_t
  }
  
  return(Theta)
}

lambda_1 <- 0.7
lambda_2 <- 0.99

# Initialize R_0 
R_0 <- diag(0.1, 2, 2)

# Calculate parameter estimates for both lambdas
Theta_lambda_1 <- RLS_estimate_with_forgetting(X, y, R_0, lambda_1)
Theta_lambda_2 <- RLS_estimate_with_forgetting(X, y, R_0, lambda_2)

# Plot the estimates
par(mfrow=c(1,2))  

# Plot for θ1
plot(1:n, Theta_lambda_1[, 1], type="l", col="blue", ylim=range(c(Theta_lambda_1[,1], Theta_lambda_2[,1])), xlab="t", ylab="θ1", main="Estimates of θ1")
lines(1:n, Theta_lambda_2[, 1], col="red")
legend("topright", legend=c("λ=0.7", "λ=0.99"), col=c("blue", "red"), lty=1)

# Plot for θ2
plot(1:n, Theta_lambda_1[, 2], type="l", col="blue", ylim=range(c(Theta_lambda_1[,2], Theta_lambda_2[,2])), xlab="t", ylab="θ2", main="Estimates of θ2")
lines(1:n, Theta_lambda_2[, 2], col="red")
legend("topright", legend=c("λ=0.7", "λ=0.99"), col=c("blue", "red"), lty=1)

# The estimates of θ1 and θ2 are more stable when using a higher value of λ.
# This is because a higher value of λ gives more weight to past observations,
# making the model less sensitive to recent fluctuations in the data.
# As a result, the estimates are smoother and more stable over time.

# On the other hand, a lower value of λ (e.g., 0.7) gives more weight to the most recent observations,
# making the model more sensitive to changes in recent data and causing the estimates to fluctuate more.
# This can lead to less stable estimates in the short term, as the model reacts strongly to recent changes.


# 4.5
# Function to compute one-step ahead predictions
one_step_forecast <- function(X, Theta) {
  n <- nrow(X)
  y_hat <- numeric(n)
  
  for (t in 5:n) {
    y_hat[t] <- sum(X[t, ] * Theta[t - 1, ])
  }
  return(y_hat)
}

# Predictions for both lambdas
y_hat_lambda_1 <- one_step_forecast(X, Theta_lambda_1)
y_hat_lambda_2 <- one_step_forecast(X, Theta_lambda_2)

# Residuals for both lambdas
residuals_lambda_1 <- y_hat_lambda_1 - y
residuals_lambda_2 <- y_hat_lambda_2 - y

# Removing the burn-in period (those are thr first 4 points)
t_values <- 5:n
residuals_lambda_1 <- residuals_lambda_1[t_values]
residuals_lambda_2 <- residuals_lambda_2[t_values]

par(mfrow=c(1,1))
plot(t_values, residuals_lambda_1, type="l", col="blue", ylim=range(c(residuals_lambda_1, residuals_lambda_2)), xlab="t", ylab="Residuals", main="One-step ahead residuals")
lines(t_values, residuals_lambda_2, col="red")
legend("topright", legend=c("λ=0.7", "λ=0.99"), col=c("blue", "red"), lty=1)


# 4.6
# Function to compute k-step-ahead predictions
k_step_forecast <- function(X, Theta, k) {
  n <- nrow(X)
  y_hat_k <- rep(NA, n)
  
  for (t in 1:(n - k)) {
    y_hat_k[t + k] <- sum(X[t, ] * Theta[t, ])
  }
  
  return(y_hat_k)
}

# Function to compute k-step residuals
compute_k_step_residuals <- function(y, y_hat_k, k) {
  valid_indices <- (k + 1):length(y)
  return(y_hat_k[valid_indices] - y[valid_indices])
}

# Function to optimize lambda across horizons k = 1,...,12
optimize_lambda <- function(X, y, lambda_seq, k_max) {
  rmse_vals <- matrix(NA, nrow = length(lambda_seq), ncol = k_max)
  colnames(rmse_vals) <- paste0("k=", 1:k_max)
  
  for (i in seq_along(lambda_seq)) {
    lambda <- lambda_seq[i]
    Theta <- RLS_estimate_with_forgetting(X, y, diag(0.1, 2, 2), lambda)
    
    for (k in 1:k_max) {
      y_hat_k <- k_step_forecast(X, Theta, k)
      residuals_k <- compute_k_step_residuals(y, y_hat_k, k)
      rmse_vals[i, k] <- sqrt(mean(residuals_k^2, na.rm = TRUE))
    }
  }
  
  return(rmse_vals)
}

# Runing the optimization over lambda
lambda_seq <- seq(0.5, 0.99, by = 0.01)
k_max <- 12
rmse_results <- optimize_lambda(X, y, lambda_seq, k_max)

# Dataframe for plotting
df <- data.frame(lambda = rep(lambda_seq, k_max), 
                 k = rep(1:k_max, each = length(lambda_seq)), 
                 RMSE = as.vector(rmse_results))

library(ggplot2)
ggplot(df, aes(x = lambda, y = RMSE, color = factor(k))) +
  geom_line() +
  labs(title = "RMSE vs Lambda", x = "Lambda", y = "RMSE", color = "Horizon k") +
  theme_minimal()




# 4.7


