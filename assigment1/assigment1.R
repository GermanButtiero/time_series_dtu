Dtrain
Dtest
#install.packages("dplyr")  # Install (only needed once)
library(dplyr)  # Load the package
library(tidyverse)

plot(Dtrain$year, Dtrain$total)

# Fit a linear model
fit <- lm(total ~ year, data = Dtrain)
summary(fit)

# abline(fit, col = "red")

# Extract the parameter estimates and their standard errors
theta_hat <- coef(fit)
theta_hat
theta_hat_se <- summary(fit)$coef[, 2]
theta_hat_se

# Add the estimated mean to the plot
lines(Dtrain$year, predict(fit), col = "red")

# Make a forecast with intervals
Dtest$forecast <- predict(fit, newdata = Dtest, interval = "prediction", level=0.95)
Dtest

ggplot(Dtrain, aes(x = time, y = total)) +
  geom_point(col="red") +
  geom_line(aes(y=predict(fit)), col="red", size=.5) +
  geom_line(data=Dtest, aes(y=forecast), col="blue", size=.5) +

