
# filter + store states
kf_filter_states <- function(par, df) {
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
#   x_est  <- matrix(Y[1,], n, 1)
  x_est <- X0
  P_est  <- diag(1e1, n)
  
  x_storage <- matrix(NA, nrow=Tn, ncol=2)
  
  for(t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t,], 3, 1)
    P_pred <- A %*% P_est %*% t(A) + Q
    
    # innovation
    y_pred <- C %*% x_pred
    S_t    <- C %*% P_pred %*% t(C) + R
    innov  <- Y[t,] - y_pred
    
    K_t   <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- (diag(n) - K_t %*% C) %*% P_pred
    
    x_storage[t, ] <- as.vector(x_est)
  }
  
  colnames(x_storage) <- c("x1", "x2")
  x_storage
}


# Load data
df <- read.csv("assignment4/transformer_data.csv")

# df$Ta_s <- scale(df$Ta)
# df$S_s  <- scale(df$S)
# df$I_s  <- scale(df$I)
# df$Y_s <- scale(df$Y)


# Estimated parameter values from 2.3
est_par <- c(
  # A: 2x2 state transition matrix
   0.8,  0.08,   # A11, A12
   0.003,  0.3,   # A21, A22
  
  # B: 2x3 input matrix (inputs are Ta, S, I)
   0.1, -0.1, 0.004,   # B11, B12, B13 (effect of Ta, S, I on state 1)
   -0.3, 0.2, 0.03,   # B21, B22, B23 (effect on state 2)
  
  # Q (process noise covariance matrix, lower-triangular Cholesky)
   0.6,   # L11
   0.0007,  # L21 
   1,   # L22
  
  # R: observation noise variance (scalar, must be > 0)
    0.001,    # R

    # X0: initial state (2x1)
    10, 10   # X01, X02
)

# est_par <- c(
#   # A: 2x2 state transition matrix
#    1.0,  -0.1,   # A11, A12
#    1.0,  -0.05,   # A21, A22
  
#   # B: 2x3 input matrix (inputs are Ta, S, I)
#    1.5, -1.0, -0.02,   # B11, B12, B13 (effect of Ta, S, I on state 1)
#    -0.4, 0.5, 0.5,   # B21, B22, B23 (effect on state 2)
  
#   # Q (process noise covariance matrix, lower-triangular Cholesky)
#    -0.005,   # L11
#    0.4,  # L21 
#    0.4,   # L22
  
#   # R: observation noise variance (scalar, must be > 0)
#     0.2,    # R

#     # X0: initial state (2x1)
#     10, 10   # X01, X02
# )

# est_par <- c(
#   # A: 2x2 state transition matrix
#    0.07,  -0.3,   # A11, A12
#    -0.1,  0.7,   # A21, A22
  
#   # B: 2x3 input matrix (inputs are Ta, S, I)
#    0.04, -0.5, 0.03,   # B11, B12, B13 (effect of Ta, S, I on state 1)
#    -0.3, 0.06, -1.5,   # B21, B22, B23 (effect on state 2)
  
#   # Q (process noise covariance matrix, lower-triangular Cholesky)
#    0.03,   # L11
#    -0.06,  # L21 
#    0.8,   # L22
  
#   # R: observation noise variance (scalar, must be > 0)
#     0.001,    # R

#     # X0: initial state (2x1)
#     10, 10   # X01, X02
# )


# Reconstructed states
states <- kf_filter_states(est_par, df)

time <- 1:nrow(df)

# Plot both states in one figure and save as PNG
png("assignment4/plots/reconstructed_states.png", width=800, height=600, res = 150, pointsize = 9)
par(mfrow = c(1,1), mar = c(4,4,2,1))

plot(time, states[,"x1"], type="l", lwd=2, col="blue",
     ylim = range(states),
     xlab="Time", ylab="State value",
     main="Reconstructed States")
lines(time, states[,"x2"], lwd=2, col="red")
legend("topright", legend=c("State 1 (x1)","State 2 (x2)"), col=c("blue","red"), lwd=2)
dev.off()


# Plot the two states + three inputs in 5 subplots and save as PNG
png("assignment4/plots/states_and_inputs.png", width=800, height=1200, res = 150, pointsize = 9)
par(mfrow = c(5,1), mar = c(2,4,2,1))

# State 1
plot(time, states[,"x1"], type="l", lwd=1.5, col="#0000FF",
     xlab="", ylab="x1", main="State 1")
# State 2
plot(time, states[,"x2"], type="l", lwd=1.5, col="red",
     xlab="", ylab="x2", main="State 2")
# Input Ta_s
plot(time, df$Ta, type="l", lwd=1.5,
     xlab="", ylab="Ta", main="Input: Ta")
# Input S_s
plot(time, df$S, type="l", lwd=1.5,
     xlab="", ylab="S", main="Input: S")
# Input I_s
plot(time, df$I, type="l", lwd=1.5,
     xlab="Time", ylab="I", main="Input: I")

dev.off()