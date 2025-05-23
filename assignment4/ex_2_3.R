
# log-likelihood only (for optim)
kf_logLik_dt <- function(par, df, return_residuals = FALSE) {
  # unpack
  A <- matrix(par[1:4], 2, 2)
  B <- matrix(par[5:10], 2, 3)
  L <- matrix(c(par[11], 0, par[12], par[13]), 2, 2)
  Q <- L %*% t(L)
  C <- matrix(c(1, 0), nrow=1)
  R <- matrix(par[14]^2, 1, 1)
  X0 <- matrix(par[15:16], 2, 1)
  
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

# Initial parameter values
start_par <- c(
  # A: 2x2 state transition matrix
   0.9,  0.0,   # A11, A12
   0.0,  0.9,   # A21, A22
  
  # B: 2x3 input matrix (inputs are Ta, S, I)
   0.01, 0.01, 0.01,   # B11, B12, B13 (effect of Ta, S, I on state 1)
   0.01, 0.01, 0.01,   # B21, B22, B23 (effect on state 2)
  
  # Q (process noise covariance matrix, lower-triangular Cholesky)
   1,   # L11
   0.01,  # L21 
   1,   # L22
  
  # R: observation noise variance (scalar, must be > 0)
    1,    # R

  # initial X0
    20.0, 20.0
)

# Lower and upper bounds for parameters
lower <- rep(-10, length(start_par))
upper <- rep(10, length(start_par))
lower[11] <- lower[13] <- 1e-3   # Q diagonal elements must be >0
lower[14] <- 1e-3                # sigma2 > 0


# Run the optimization
result_2d <- estimate_dt(start_par, df, lower, upper)
params_2d <- result_2d$par

residuals_2d <- kf_logLik_dt(params_2d, df, return_residuals = TRUE)

pdf("plots/exer23.pdf", width = 10, height = 6)
par(mfrow = c(2, 2))

plot(residuals_2d, type = 'l', main = 'Residuals', ylab = "Residuals", xlab = "Time")
acf(residuals_2d, main = 'ACF of Residuals')
pacf(residuals_2d, main = 'PACF of Residuals')
qqnorm(residuals_2d); qqline(residuals_2d)

logLik_2d <- -kf_logLik_dt(params_2d, df)
n <- nrow(df)
k <- length(params_2d)
AIC_2d <- -2 * logLik_2d + 2 * k
BIC_2d <- -2 * logLik_2d + log(n) * k
cat("AIC (2D):", AIC_2d, "\nBIC (2D):", BIC_2d, "\n")

dev.off()
