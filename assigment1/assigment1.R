#install.packages("dplyr")  # Install (only needed once)
library(dplyr)  # Load the package
library(tidyverse)

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