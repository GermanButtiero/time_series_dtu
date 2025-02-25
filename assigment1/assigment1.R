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
  geom_line(aes(y = predict(fit)), col = "red", size = 0.8) +
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

# 4.1
