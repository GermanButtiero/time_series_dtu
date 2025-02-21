
##############################################
### Examples for Time Series Aanalysis 02417, Lecture 3 (2025)
# 1 - Linear regression with matrix notation 
# 2 - Weighted Least Squares
# 3 - From Global to Local model (OLS with fewer datapoints)
# 4 - WLS with "lambda weights" (exponential weights)
# 5 - RLS with forgetting
##############################################


# load some packages:
library(fpp2)
library(dplyr)
library(tidyverse)


# load data:
# The "austa" data = Total international visitors to Australia (in millions) 
# from 1980 to 2015. Yearly data (ie one datapoint per year). In total 36 datapoints.

autoplot(austa)
# austa is a "time-series" object, but we will convert to a "dataframe" and call it austa_data
austa_data <- data.frame(int_visit=as.matrix(austa), year=as.matrix(time(austa)))
rownames(austa_data) <- NULL
head(austa_data)

# Now we split the data into test and train:
austa_train <- austa_data %>% filter(., year <= 2005)
austa_test <- austa_data %>% filter(., year > 2005)


### Plot train and test data
ggplot(austa_train, aes(x=year, y=int_visit)) +
  geom_point(col="red") + 
  geom_point(data=austa_test, col="black") +
  xlim(1980, 2020) + ylim(0, 8)


#################################################
# 1- Linear regression with matrix notation
#################################################

# we will fit a linear model: Y_i = beta_0 + beta_1 * time_i + epsilon_i
# so we have two parameters: beta_0 and beta_1
p <- 2

# we also save the number of observations (26 obs. in the training data):
n <- length(austa_train$year)

# X is the "design matrix"
X <- cbind(1, austa_train$year)
print(X)

# y is vector with observations:
y <- cbind(austa_train$int_visit)
print(y)

# to estimate parameters we solve the "normal equations":
theta_OLS <- solve(t(X)%*%X)%*%t(X)%*%y
print(theta_OLS)
# these are the parameter estimates!
theta_0 <- theta_OLS[1]
theta_1 <- theta_OLS[2]

# now we compute y_hat values (so thet we can plot the regression line) 
yhat_ols <- X%*%theta_OLS


# plot:
ggplot(austa_train, aes(x=year, y=int_visit)) +
  geom_point(col="red") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5) +
  xlim(1980, 2020) + ylim(0, 8)


# we will now calculate the standard errors on the parameters beta_0 and beta_1:

# first compute residuals:
e_ols <- y - yhat_ols

# calculate sum of squared residuals:
RSS_ols <- t(e_ols)%*%e_ols

# calculate sigma^2:
sigma2_ols <- as.numeric(RSS_ols/(n - p))

# calculate variance-covariance matrix of _parameters_:
V_ols <- sigma2_ols * solve(t(X) %*% X)
print(V_ols)

# the variances of the parameters are the values in the diagonal:
diag(V_ols)
# and the standard errors are given by:
sqrt(diag(V_ols))

se_theta_0 <- (sqrt(diag(V_ols)))[1] # standard error of the intercept-parameter
se_theta_1 <- (sqrt(diag(V_ols)))[2] # standard error of the slope-parameter

# now we have both point estimates and standard errors:
# intercept:
theta_0
se_theta_0
# slope:
theta_1
se_theta_1

# Predictions for future values ("forecast")
# now we use the model for predictions on future timepoints
# we use the timepoints from the testdata:
Xtest <- cbind(1, austa_test$year)
print(Xtest)

# compute predictions (we compute all 10 predictions at once):
y_pred <- Xtest%*%theta_OLS
print(y_pred)

# compute prediction variance-covariance matrix:
Vmatrix_pred <- sigma2_ols*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
# the variances of individual predictions are in the diagonal of the matrix above
print(diag(Vmatrix_pred))

# compute "prediction intervals" 
y_pred_lwr <- y_pred - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
y_pred_upr <- y_pred + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))


# plot forecast:
ggplot(austa_train, aes(x=year, y=int_visit)) +
  geom_point(col="red") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5) +
  geom_point(data=austa_test, aes(x=year,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=austa_test, aes(x=year,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.2, fill="red") +
  xlim(1980, 2020) + ylim(0, 8)
austa_test
# plot WITH true test data:
ggplot(austa_train, aes(x=year, y=int_visit)) +
  geom_point(col="red") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5) +
  geom_point(data=austa_test, aes(x=year,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=austa_test, aes(x=year,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.2, fill="red") +
  geom_point(data=austa_test, aes(x=year,y=int_visit), col="black") +
  xlim(1980, 2020) + ylim(0, 8)


# (back to slides)

# qq plot of residuals:
qqnorm(e_ols)
qqline(e_ols)

# plot residuals versus x (year):
ggplot(austa_train, aes(x=year)) +
  geom_point(aes(y=e_ols), col="blue") +
  geom_line(aes(y=e_ols), col="blue") + 
  ylim(-1,1)

# plot some white noise:
set.seed(876573)
white_noise = rnorm(n=26, mean = 0, sd = sqrt(sigma2_ols))
qqnorm(white_noise)
qqline(white_noise)

ggplot(austa_train, aes(x=year)) +
  geom_point(aes(y=white_noise), col="blue") +
  geom_line(aes(y=white_noise), col="blue") + 
  ylim(-1, 1)

# (back to slides)


#################################################
# 2 - Weighted Least Squares
#################################################

# set up SIGMA matrix:
SIGMA <- diag(n)
rho <- 0.6 
SIGMA[1:5,1:5]
for (row in 1:n) {
  for (col in 1:n) {
    if(row==col) {
      SIGMA[row,col] <- 1
    } else {
      SIGMA[row,col] <- rho^abs(row-col)
    }
  }
}
# print upper corner to check:
SIGMA[1:5,1:5] #looks fine :-)


# estimate parameters with WLS
theta_WLS <- solve(t(X)%*%solve(SIGMA)%*%X)%*%(t(X)%*%solve(SIGMA)%*%y)
print(theta_WLS)

yhat_wls <- X%*%theta_WLS
e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%solve(SIGMA)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - p))

V_wls <- sigma2_wls *  solve(t(X)%*%solve(SIGMA)%*%X)
print(sqrt(diag(V_wls))) 

# predictions with WLS estimates:
y_pred_wls <- Xtest%*%theta_WLS
Vmatrix_pred <- sigma2_wls * (1 + (Xtest %*% solve(t(X)%*%solve(SIGMA)%*%X)) %*% t(Xtest) )
y_pred_lwr_wls <- y_pred_wls - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
y_pred_upr_wls <- y_pred_wls + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))

# plot WLS (blue) together with OLS (red) and true test data (black):
ggplot(austa_train, aes(x=year, y=int_visit)) +
  geom_point(col="red") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5) +
  geom_point(data=austa_test, aes(x=year,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=austa_test, aes(x=year,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.2, fill="red") +
  geom_point(data=austa_test, aes(x=year,y=int_visit), col="black") +
  xlim(1980, 2020) + ylim(0, 8) + 
  geom_line(aes(y=yhat_wls), col="blue", size=.5) +
  geom_point(data=austa_test, aes(x=year,y=y_pred_wls), col="blue", size=.5) +
  geom_ribbon(data=austa_test, aes(x=year,ymin=y_pred_lwr_wls, ymax=y_pred_upr_wls), inherit.aes=FALSE, alpha=0.2, fill="blue") 

  
# Was our guess of rho a good guess?  
# lets calculate the correlation between lag-1 residuals:
cor(e_wls[1:(n-1)], e_wls[2:n], method="pearson")
# maybe we could try with rho = 0.7 

# (back to slides)

#################################################
# 3 - From Global to Local model (OLS with fewer datapoints)
#################################################

# Fit data with less datapoints:
fit_austa <- function(start, stop) {
  n <- stop - start + 1
  p <- 2
  theta_OLS <- solve(t(X[start:stop,])%*%X[start:stop,])%*%t(X[start:stop,])%*%y[start:stop]
  
  yhat <- X%*%theta_OLS
  Xpred <- cbind(1, seq(1980+stop, 1980+stop+9, 1))
  ypred <- Xpred%*%theta_OLS
  pred <- data.frame(Xpred, ypred)
  
  RSS_ols <- sum((yhat[start:stop] - y[start:stop])^2)
  sigma2_ols <- as.numeric(RSS_ols/(n - p))
  Vmatrix_pred <- sigma2_ols*(1+(Xpred%*%solve(t(X[start:stop,])%*%X[start:stop,]))%*%t(Xpred))
  pred$ypred_lwr <- pred$ypred - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
  pred$ypred_upr <- pred$ypred + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
  
  myplot <- ggplot(austa_data, aes(x=year, y=int_visit)) +
    geom_point() + 
    geom_point(data=austa_data[start:stop,], col="red") + 
    geom_line(data=austa_data[1:(stop),], aes(y=yhat), col="red", size=.5) + 
    geom_point(data=pred, aes(x=Xpred[,2], y=ypred), col="red", size=.5) + 
    geom_ribbon(data=pred, aes(x=Xpred[,2],ymin=ypred_lwr, ymax=ypred_upr), inherit.aes=FALSE, alpha=0.2, fill="red") +
    xlim(1980, 2020) + ylim(0, 8)
  
  print(myplot)
}

fit_austa(1,26)
fit_austa(17,26)
fit_austa(21,26)
fit_austa(25,26)




#################################################
# 4 - WLS with "lambda weights" (exponential weights)
#################################################

lambda = 0.6
weights <- lambda^((n-1):0)
# plot the weights:
barplot(weights, names=1:26)


SIGMA <- diag(n)
diag(SIGMA) <- 1/weights
# print lower right corner to check:
print(SIGMA[20:26,20:26]) # looks fine :-)


# estimate parameters with WLS
theta_WLS <- solve(t(X)%*%solve(SIGMA)%*%X)%*%(t(X)%*%solve(SIGMA)%*%y)
print(theta_WLS)
yhat_wls <- X%*%theta_WLS

ggplot(austa_train, aes(x=year, y=int_visit)) +
  geom_point(col="black") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5, linetype=2) +
  geom_line(aes(y=yhat_wls), col="blue", size=.5)

# (back to slides)

# predictions from the WLS model:
y_pred_wls <- Xtest%*%theta_WLS

# prediction intervals:
e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%solve(SIGMA)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - p))
Vmatrix_pred <- sigma2_wls * (1 + (Xtest %*% solve(t(X)%*%solve(SIGMA)%*%X)) %*% t(Xtest) )
y_pred_lwr_wls <- y_pred_wls - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
y_pred_upr_wls <- y_pred_wls + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))

# plot WLS (blue) together with OLS (red) and true test data (black):
ggplot(austa_train, aes(x=year, y=int_visit)) +
  geom_point(col="red") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5, linetype=2) +
  geom_point(data=austa_test, aes(x=year,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=austa_test, aes(x=year,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.1, fill="red") +
  geom_point(data=austa_test, aes(x=year,y=int_visit), col="black") +
  geom_line(aes(y=yhat_wls), col="blue", size=.5) +
  geom_point(data=austa_test, aes(x=year,y=y_pred_wls), col="blue", size=.5) +
  geom_ribbon(data=austa_test, aes(x=year,ymin=y_pred_lwr_wls, ymax=y_pred_upr_wls), inherit.aes=FALSE, alpha=0.2, fill="blue") +
  coord_cartesian(ylim = c(0, 8), xlim = c(1980,2020)) 

# these prediction intervals are very (too) narrow!

# better prediction intervals:
Tmemory <- sum(1/diag(SIGMA))
print(Tmemory)

sigma2_wls <- as.numeric(RSS_wls/(Tmemory - p))
Vmatrix_pred <- sigma2_wls * (1 + (Xtest %*% solve(t(X)%*%solve(SIGMA)%*%X)) %*% t(Xtest) )
y_pred_lwr_wls <- y_pred_wls - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))
y_pred_upr_wls <- y_pred_wls + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred))

# plot WLS (blue) together with OLS (red) and true test data (black):
ggplot(austa_train, aes(x=year, y=int_visit)) +
  geom_point(col="red") + 
  geom_line(aes(y=yhat_ols), col="red", size=.5, linetype=2) +
  geom_point(data=austa_test, aes(x=year,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=austa_test, aes(x=year,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.1, fill="red") +
  geom_point(data=austa_test, aes(x=year,y=int_visit), col="black") +
  geom_line(aes(y=yhat_wls), col="blue", size=.5) +
  geom_point(data=austa_test, aes(x=year,y=y_pred_wls), col="blue", size=.5) +
  geom_ribbon(data=austa_test, aes(x=year,ymin=y_pred_lwr_wls, ymax=y_pred_upr_wls), inherit.aes=FALSE, alpha=0.2, fill="blue") +
  coord_cartesian(ylim = c(0, 8), xlim = c(1980,2020)) 




#################################################
# 5 - RLS with forgetting
#################################################

# sanity check:
# using RLS notation for computing OLS with all 26 observations:
R26 <- t(X)%*%X
print(R26)

h26 <- t(X)%*%y
print(h26)

RLS26 <- solve(R26)%*%h26
print(RLS26)
# this corresponds exactly to the OLS estimates (as it should)

# (back to slides)


###########################
# Now iterate through data:
# For each step:
# - calculate R and theta
# - calculate one-step prediction

lambda <- 0.6

# we will use the entire dataset (not only the training data from before):
X <- cbind(1, austa_data$year)
y <- cbind(austa_data$int_visit)

n <- length(X[,1])

# initialise containers for parameter estimates (Theta) and one step predictions:
Theta <- matrix(NA, nrow=n, ncol=p)
OneStepPred <- matrix(NA, nrow=n)

# 1 # very first step:
x1 <- X[1,]

R_1 <- x1%*%t(x1) # R is a pxp matrix
h_1 <- x1*y[1]    # h is a px1 vector (but R prints it in a row..)

# to estimate theta we need to invert R:
solve(R_1)
# in this very first step R cannot be inverted - too soon to estimate parameters!
# (we cannot estimate p parameters drom only one datapoint)


# 2 # second step - first time to estimate parameters and make prediction
x2 <- X[2,]
R_2 <- lambda*R_1 + x2 %*% t(x2)
h_2 <- lambda*h_1 + x2 * y[2]

solve(R_2)
# R is now invertible (we can estimate p parameters from p observations)

# we estimate theta (for the first time - so not yet using "update" formula):
Theta[2,] <- solve(R_2) %*% h_2

# we predict one step ahead:
OneStepPred[2+1] <- X[2+1,]%*%Theta[2,]


# 3 # third step - first time to use update formula
x3 <- X[3,]
R_3 <- lambda*R_2 + x3 %*% t(x3)
Theta[3,] <- Theta[2,] + solve(R_3) %*% x3 %*% (y[3] - t(x3) %*% Theta[2,])

# we predict one step ahead:
OneStepPred[3+1] <- X[3+1,]%*%Theta[3,]

# next many steps # - update and predict

R <- R_3

for(i in 4:n){
  x <- X[i, ]
  # Update
  R <- lambda*R + x %*% t(x)
  Theta[i, ] <- Theta[i-1, ] + solve(R) %*% x %*% (y[i] - t(x) %*% Theta[i-1, ])
}

# predict
for(i in 4:n-1){
  OneStepPred[i+1] <- X[i+1, ]%*%Theta[i, ]
}


# Plot estimate of intercept:
plot(Theta[,1])

# Plot estimate of slope:
plot(Theta[,2])

# Plot one step predictions:
ggplot(austa_data, aes(x=year, y=int_visit)) +
  geom_point(col="black") +
  geom_point(aes(y=OneStepPred), col="blue", size=1) + 
  geom_line(aes(y=OneStepPred), col="blue") + 
  coord_cartesian(ylim = c(0, 8), xlim = c(1980,2020)) 



