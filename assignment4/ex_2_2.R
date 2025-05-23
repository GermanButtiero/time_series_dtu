
# log-likelihood only (for optim)
kf_logLik_dt <- function(par, df, return_residuals = FALSE) {
  # unpack
  A <- matrix(par[1], 1, 1)
  B <- matrix(par[2:4], 1, 3)
  L <- matrix(par[5], 1, 1) 
  Q <- L %*% t(L)
  C <- matrix(1, 1, 1)
  R <- matrix(par[6]^2, 1, 1)
  X0 <- matrix(par[7], 1, 1)
  
  # data
  Y <- as.matrix(df[,"Y",drop=FALSE])
  U <- as.matrix(df[,c("Ta","S","I")])
  Tn <- nrow(df)
  
  # initialize
  n      <- nrow(A)
  x_est  <- matrix(Y[1,], n, 1)
  x_est <- X0
  P_est  <- diag(1e1, n)
  logLik <- 0

  residuals <- numeric(Tn)
  
  for(t in 1:Tn) {
    # predict
    x_pred <- A %*% x_est + B %*% t(U[t, , drop = FALSE])
    P_pred <- A %*% P_est %*% t(A) + Q
    
    # innovation
    y_pred <- C %*% x_pred
    S_t    <- C %*% P_pred %*% t(C) + R
    innov  <- Y[t,] - y_pred

    residuals[t] <- innov  # store residual
    
    # accumulate log-lik
    logLik <- logLik - 0.5 * (log(2 * pi) + log(det(S_t)) + t(innov) %*% solve(S_t, innov))
    
    # update
    K_t   <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- (diag(n) - K_t %*% C) %*% P_pred
  }
  
  if (return_residuals) {
    return(residuals)
  } else {
    return(as.numeric(logLik))
  }

  as.numeric(logLik)
}

# Optimizer wrapper
estimate_dt <- function(start_par, df, lower=NULL, upper=NULL) {
  negLL <- function(par){ -kf_logLik_dt(par, df) }
  optim(
    par    = start_par, fn = negLL,
    method = "L-BFGS-B",
    lower  = lower, upper = upper,
    control= list(maxit=1000, trace=1)
  )
}




### Load data
df <- read.csv("assignment4/transformer_data.csv")
df

# Initial parameter values
start_par <- c(
  0.9,         # A
  0.01, 
  0.01, 
  0.01,  # B1, B2, B3
  1.0,         # L (Cholesky of Q)
  1.0,          # R (observation noise variance)
  20.0          # initial X0
)

# Lower and upper bounds for parameters
lower <- c(-1.5, -1, -1, -1, 1e-3, 1e-3, 10)
upper <- c(1.5, 1, 1, 1, 10, 10, 30)

# Run the optimization
result <- estimate_dt(start_par, df, lower, upper)
params <- result$par

# Residuals
residuals <- kf_logLik_dt(params, df, return_residuals = TRUE)

# Plots
par(mfrow=c(2,2))
plot(residuals, type='l', main='Residuals')
acf(residuals, main='ACF of Residuals')
pacf(residuals, main='PACF of Residuals')
qqnorm(residuals); qqline(residuals)

# AIC and BIC
logLik <- -kf_logLik_dt(params, df)
n <- nrow(df)
k <- length(params)
AIC <- -2 * logLik + 2 * k
BIC <- -2 * logLik + log(n) * k
cat("AIC:", AIC, "\nBIC:", BIC, "\n")

# Residual diagnostics
pdf("images/exer22.pdf", width = 10, height = 6)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Residuals plot
plot(residuals, type = 'l', main = 'Residuals', ylab = "Residuals", xlab = "Time")

# ACF plot
acf(residuals, main = 'Autocorrelation Function (ACF)')

# PACF plot
pacf(residuals, main = 'Partial Autocorrelation Function (PACF)')

# QQ plot
qqnorm(residuals, main = 'Normal Q-Q Plot')
qqline(residuals)

dev.off()
