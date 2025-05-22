
kf_logLik_dt <- function(par, df) {
  # par: vector of parameters
  # df: data frame with observations and inputs as columns (Y, Ta, S, I)
  # par: Could be on the form c(A11, A12, A21, A22, B11, B12, B21, B22, Q11, Q12, Q22)

  A   <- matrix(par[1:4], 2,2, byrow=TRUE) # transition matrix
  B   <- matrix(par[5:10], 2,3, byrow=TRUE) # input matrix

  L <- matrix(c(par[11], 0, par[12], par[13]), 2,2, byrow=TRUE) # lower-triangle of system covariance matrix

  Q <- L %*% t(L) # THAT IS!!! The system covariance matrix is given by Qlt %*% t(Qlt) (and is this symmetric positive definite)

  C   <- matrix(c(1, 0), 1,2) # observation matrix
  R <- matrix(par[14], 1,1) # observation noise covariance matrix

  X0  <- matrix(df$Y[1], 2,1) # initial state

  # Variables
  obs_cols <- c("Y") # observation column names
  input_cols <- c("Ta","S","I") # input column names

  # pull out data
  Y  <- as.matrix(df[, obs_cols])     # m×T
  U  <- as.matrix(df[, input_cols])   # p×T

  # init
  Tn     <- nrow(df)
  x_est  <- matrix(Y[1,], 2, 1)            # start state from first obs
  # x_est  <- X0 
  P_est  <- diag(1e1, 2)                   # X0 prior covariance
  logLik <- 0

  for (t in 1:Tn) {
    # prediction step
    x_pred <- A %*% x_est + B %*% matrix(U[t,], 3,1) # write the prediction step
    P_pred <- A %*% P_est %*% t(A) + Q # write the prediction step (Sigma_xx)

    # innovation step
    y_pred  <- C %*% x_pred # predicted observation
    S_t     <- C %*% P_pred %*% t(C) + R # predicted observation covariance (Sigma_yy)
    innov   <- matrix(Y[t,],1,1) - y_pred # innovation (one-step prediction error)

    # log-likelihood contribution
    logLik <- logLik - 0.5*(sum(log(2*pi*S_t)) + t(innov) %*% solve(S_t, innov))

    # update step
    K_t    <- P_pred %*% t(C) %*% solve(S_t) # Kalman gain
    x_est  <- x_pred + K_t %*% innov # reconstructed state
    P_est  <- (diag(2) - K_t %*% C) %*% P_pred # reconstructed covariance
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
  # A: 2x2 state transition matrix
   0.9,  0.1,   # A11, A12
   0.0,  0.9,   # A21, A22
  
  # B: 2x3 input matrix (inputs are Ta, S, I)
   0.05, 0.01, 0.1,   # B11, B12, B13 (effect of Ta, S, I on state 1)
   0.02, 0.01, 0.1,   # B21, B22, B23 (effect on state 2)
  
  # Q (process noise covariance matrix, lower-triangular Cholesky)
   0.1,   # L11
   0.01,  # L21
   0.1,   # L22
  
  # R: observation noise variance (scalar, must be > 0)
   0.5    # R
)

# Lower and upper bounds for parameters
lower <- c(rep(-1, 4),   # A entries
           rep(-1, 6),   # B entries
           1e-4, -1, 1e-4, # L11, L21, L22 (L11, L22 > 0 for Cholesky)
           1e-4)          # R > 0

upper <- c(rep(1.5, 4),   # A
           rep(1, 6),     # B
           10,  10, 10,   # Q
           10)            # R

# Run the optimization
estimate_dt(start_par, df, lower = lower, upper = upper)
