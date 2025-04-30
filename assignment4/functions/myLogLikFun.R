source("assignment4/functions/kalmanfilter.R")

myLogLikFun <- function(theta, y, R, x_prior = 0, P_prior = 10) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  
  kf_result <- myKalmanFilter(y, theta, R, x_prior, P_prior) # call the Kalman filter function
  err <- kf_result$innovation       # Innovations
  S <- kf_result$innovation_var   # Innovation covariances
  
  # Compute log-likelihood contributions from each time step
  n <- length(y)
  logL_components <- -0.5 * (log(2 * pi) + log(S) + (err^2 / S))
  logL <- sum(logL_components)
  
  return(-logL)  # Return negative log-likelihood for minimization
}
