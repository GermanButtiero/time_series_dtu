Dtrain
Dtest
#install.packages("dplyr")  # Install (only needed once)
library(dplyr)  # Load the package

create_x <- function(df) {
  df$time <- as.Date(df$time, format = "%Y-%m-%d")
  df$year <- as.numeric(format(df$time, "%Y"))
  df$month <- as.numeric(format(df$time, "%m"))
  df$x <- df$year + (df$month - 1) / 12
  return(df)
}
Dtrain <- create_x(Dtrain)

plot(Dtrain$x, Dtrain$total)

# Fit a linear model
fit <- lm(total ~ x, data = Dtrain)
summary(fit)

# abline(fit, col = "red")

# Extract the parameter estimates and their standard errors
theta_hat <- coef(fit)
theta_hat
theta_hat_se <- summary(fit)$coef[, 2]
theta_hat_se

# Add the estimated mean to the plot
lines(Dtrain$x, predict(fit), col = "red")

Dtest <-create_x(Dtest)
Dtest$x

# Make a forecast with intervals
Dtest$forecast <- predict(fit, newdata = Dtest, interval = "prediction", level=0.95)
Dtest

#append the time and forecast to a new data frame
times <- Dtrain$time %>% c(Dtest$time)
times
forecasts <- Dtrain$total %>% c(Dtest$forecast[, 1])
forecasts
plot(times, forecasts)
