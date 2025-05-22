
kf_logLik_dt <- function(par, df) {
  # par: vector of parameters
  # df: data frame with observations and inputs as columns (Y, Ta, S, I)
  # par: Could be on the form c(A11, A12, A21, A22, B11, B12, B21, B22, Q11, Q12, Q22)

  A   <- par[1] # transition matrix
  B   <- matrix(par[2:4], 1, 3, byrow=TRUE) # input matrix

  L <- par[5] # lower-triangle of system covariance matrix

  Q <- L %*% t(L) # THAT IS!!! The system covariance matrix is given by Qlt %*% t(Qlt) (and is this symmetric positive definite)

  C   <- 1 # observation matrix
  R <- par[6] # observation noise covariance matrix

  # Variables
  obs_cols <- c("Y") # observation column names
  input_cols <- c("Ta","S","I") # input column names

  # pull out data
  Y  <- as.matrix(df[, obs_cols])     # m×T
  U  <- as.matrix(df[, input_cols])   # p×T

  # init
  Tn     <- nrow(df)
  x_est  <- matrix(Y[1,], 1, 1)            # start state from first obs
  P_est  <- matrix(1e1, 1, 1)                  # X0 prior covariance
  logLik <- 0

  for (t in 1:Tn) {
    # prediction step
    x_pred <- A * x_est + B %*% matrix(U[t,], 3, 1) # write the prediction step
    P_pred <- A * P_est * A + Q # write the prediction step (Sigma_xx)

    # innovation step
    y_pred  <- C * x_pred # predicted observation
    S_t     <- C * P_pred * C + R # predicted observation covariance (Sigma_yy)
    innov   <- Y[t, ] - y_pred # innovation (one-step prediction error)

    # log-likelihood contribution
    logLik <- logLik - 0.5*(sum(log(2*pi*S_t)) + t(innov) %*% solve(S_t, innov))

    # update step
    K_t    <- P_pred * C / S_t # Kalman gain
    x_est  <- x_pred + K_t * innov # reconstructed state
    P_est  <- (1 - K_t * C) * P_pred # reconstructed covariance
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
  0.05, 0.01, 0.1,  # B1, B2, B3
  0.1,         # L (Cholesky of Q)
  0.5          # R (observation noise variance)
)

# Lower and upper bounds for parameters
lower <- c(
  -1.5,         # A
  rep(-1, 3),   # B
  1e-4,         # L > 0
  1e-4          # R > 0
)

upper <- c(
  1.5,          # A
  rep(1, 3),    # B
  10,           # L
  10            # R
)

# Run the optimization
estimate_dt(start_par, df, lower = lower, upper = upper)
