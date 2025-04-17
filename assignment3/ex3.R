data <- read.csv("assignment3/box_data_60min.csv")

#3.1
par(mfrow = c(3, 1))
plot(data$thour, data$Ph, type = "l", col = "blue", xlab = "Time (hours)", ylab = "Ph", main = "The heat from electrical heaters vs Time")
plot(data$thour, data$Tdelta, type = "l", col = "red", xlab = "Time (hours)", ylab = "Tdelta", main ="The difference between the internal and external temperature vs Time")
plot(data$thour, data$Gv, type = "l", col = "red", xlab = "Time (hours)", ylab = "Gv", main = "The vertical solar radiation onto the box side with a window  vs Time")
#By looking at the plots it seems that there are some dependencies between the Ph and Gv variables.
#In particular, the Ph variable seems to be negatively affected by the Gv variable.

#3.2
end_of_training <- as.POSIXct("2013-02-06 00:00:00")
#split the data into training and test sets
train <- data[data$tdate <= end_of_training, ]
test <- data[data$tdate > end_of_training, ]

nrow(train) 
nrow(test) 
#3.3
#investigate relationship between the variables with plots
# 1. Basic scatter plots to investigate relationships between Ph and other variables
par(mfrow = c(3, 1))
plot(train$Gv, train$Ph, xlab = "Vertical Solar Radiation (Gv)", ylab = "Heat from Electrical Heaters (Ph)", 
     main = "Ph vs Gv", col = "blue", pch = 16)
plot(train$Tdelta, train$Ph, xlab = "Temperature Difference (Tdelta)", ylab = "Heat from Electrical Heaters (Ph)", 
     main = "Ph vs Tdelta", col = "red", pch = 16)
plot(train$Gv, train$Tdelta, xlab = "Vertical Solar Radiation (Gv)", ylab = "Temperature Difference (Tdelta)",
     main = "Tdelta vs Gv", col = "green", pch = 16)
 
# 2. Auto-correlation of Ph to examine temporal dependencies
par(mfrow = c(1, 1))
acf(train$Ph, lag.max = 24, main = "Autocorrelation of Ph")
#it looks like the Ph variable is very dependent of the previous values
# 3. Cross-correlation between Ph and other variables
par(mfrow = c(2, 1))
ccf(train$Ph, train$Gv, main = "Cross-correlation: Ph vs Gv", , xlim = c(0, 10))
ccf(train$Ph, train$Tdelta, main = "Cross-correlation: Ph vs Tdelta", , xlim = c(0, 10))
#the ccf doesnt say much because the data is not just white noise

#3.4
formula_delta <- "Ph ~ 0 + Tdelta.l0 + Tdelta.l1 + Tdelta.l2 + Tdelta.l3 + 
Tdelta.l4 + Tdelta.l5 + Tdelta.l6 + Tdelta.l7 + Tdelta.l8 + Tdelta.l9 + Tdelta.l10"

formula_gv <- "Ph ~ 0 + Gv.l0 + Gv.l1 + Gv.l2 + Gv.l3 +
Gv.l4 + Gv.l5 + Gv.l6 + Gv.l7 + Gv.l8 + Gv.l9 + Gv.l10"

par(mfrow = c(2, 1))
model <- lm(formula_delta, data = train)
summary(model)

plot(0:10, model$coefficients, type = "h", xlab = "Lag", ylab = "Coefficient", 
     main = "Impulse Response Function: Tdelta")
model_gv <- lm(formula_gv, data = train)

summary(model_gv)
plot(0:10, model_gv$coefficients, type = "h", xlab = "Lag", ylab = "Coefficient", 
     main = "Impulse Response Function: Gv")
 
 #3.5
model_full <- lm(Ph ~ 0+ Tdelta + Gv, data = train)
summary(model_full)
#one step prediction
predictions <- predict(model_full, newdata = test)
#plot the predictions
par(mfrow = c(3, 1))
plot(test$thour, test$Ph, type = "l", col = "blue", xlab = "Time (hours)", ylab = "Ph", main = "Predictions vs Actual")
lines(test$thour, predictions, col = "red")

residuals <- test$Ph - predictions
plot(test$thour, residuals, type = "l", col = "blue", xlab = "Time (hours)", ylab = "Residuals", main = "Residuals vs Time")
acf(residuals, main = "Autocorrelation of Residuals")

#3.6
# 1. Fit a linear regression model with Tdelta and Gv as predictors
model_full <- lm(Ph ~ 0 + Ph.l1 + Tdelta + Gv, data = train)
summary(model_full)
# 2. One-step prediction
predictions <- predict(model_full, newdata = test)
# 3. Plot the predictions
par(mfrow = c(3, 1))
plot(test$thour, test$Ph, type = "l", col = "blue", xlab = "Time (hours)", ylab = "Ph", main = "Predictions vs Actual")
lines(test$thour, predictions, col = "red")
# 4. Residuals
residuals <- test$Ph - predictions
# 5. Plot the residuals
plot(test$thour, residuals, type = "l", col = "blue", xlab = "Time (hours)", ylab = "Residuals", main = "Residuals vs Time")
# 6. ACF and CCF of residuals
acf(residuals, main = "Autocorrelation of Residuals")

#3.7
model_full <- lm(Ph ~ 0+ Ph.l1 + Ph.l2 +Tdelta +Tdelta.l1 + Gv + Gv.l1, data = train)
summary(model_full)
# 2. One-step prediction
predictions <- predict(model_full, newdata = test)
# 3. Plot the predictions
par(mfrow = c(3, 1))     
plot(test$thour, test$Ph, type = "l", col = "blue", xlab = "Time (hours)", ylab = "Ph", main = "Predictions vs Actual")
lines(test$thour, predictions, col = "red")
# 4. Residuals
residuals <- test$Ph - predictions
# 5. Plot the residuals  
plot(test$thour, residuals, type = "l", col = "blue", xlab = "Time (hours)", ylab = "Residuals", main = "Residuals vs Time")
# 6. ACF and CCF of residuals
acf(residuals, main = "Autocorrelation of Residuals")

