data <- read.csv("assignment3/box_data_60min.csv")
#functions
validate <- function(fit){
  if(class(fit)[1] == "lm"){
    i <- as.integer(names(fit$residuals))
    res <- rep(NA,max(i))
    res[i] <- fit$residuals
  }else if(class(fit)[1] == "marima"){
    res <- fit$residuals[1, ]
  }else if(class(fit)[1] == "Arima"){
    res <- fit$residuals
  }
  layout(rbind(1,2:3))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7,0))
  plot(res, type="b")
  legend("topright", c("Residuals"), lty=1, col=1:2)
  acf(res, na.action=na.pass)
  title(main="ACF(residuals)", line=0.2)
  pacf(res, na.action=na.pass)
  title(main="PACF(residuals)", line=0.2)
  # Return the residuals
  invisible(res)
}

plotit <- function(predictions, residuals) {
  layout(rbind(1, c(2, 3)))
  par(mar = c(4, 4, 3, 2))  
  plot(test$thour, test$Ph, type = "l", col = "blue", 
       xlab = "Time (hours)", ylab = "Ph", 
       main = "Predictions vs Actual")
  lines(test$thour, predictions, col = "red")
legend("topright", legend = c("Actual", "Predicted"), 
         col = c("blue", "red"), lty = 1)
  acf(residuals, main = "Autocorrelation of Residuals")
  pacf(residuals, main = "Partial Autocorrelation of Residuals")
}
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

# 3.5 Fit a linear regression model with Tdelta and Gv as predictors

model_order0<- lm(Ph ~ Tdelta + Gv, data = train)
summary(model_order0)
validate(model_order0)

# The time plot of residuals shows clear patterns and cycles rather than random noise. 
# Large spikes are visible at multiple points.

# The ACF plot reveals significant autocorrelation structure, the first lag shows 
# very high correlation (approximately 0.9), followed by gradually decaying but 
# significant correlations through lag 5. This exponentially decaying pattern is 
# a classic signature of an AR process.

# The PACF plot shows a very large spike at lag 1 (around 0.7) with only minor spikes 
# afterward, suggesting an AR(1) process might be appropriate for modeling the error structure.
# 
# These diagnostics strongly indicate that a transfer function model is needed. The current 
# linear regression assumes independent errors, but the strong temporal dependence in residuals 
# means predictions would be systematically biased. A transfer function incorporating an AR(1) 
# error structure would likely provide much better forecasts by accounting for this temporal 
# correlation pattern.

#3.6
model_order1<- lm(Ph ~ Ph.l1 + Tdelta + Gv, data = train)
summary(model_order1)

validate(model_order1)
# The time plot of residuals now shows much more random behavior, though there 
# are still some notable spikes (around indices 90, 100). Overall, the residuals
# appear more stationary than in the previous model.
#
# The ACF plot is dramatically improved compared to the original model. There is
# no longer the strong autocorrelation pattern. Only the first few lags show
# minor significance, which is a substantial improvement.
#
# The PACF plot shows smaller spikes compared to the original model, with only
# a few scattered significant lags. This suggests much of the AR structure has
# been captured by including Ph.l1 in the model.
#
# Adding the lagged dependent variable (Ph.l1) has effectively addressed much of
# the temporal dependence issue seen in the original model. This is essentially
# a simplified transfer function model with an AR(1) component incorporated directly
# into the regression.

#3.7
model_order2 <- lm(Ph ~ Ph.l1 + Ph.l2 +Tdelta +Tdelta.l1 + Gv + Gv.l1, data = train)
summary(model_order2)

validate(model_order2)
# The time plot of residuals from this expanded model with multiple lags shows 
# a much more random pattern. While there are still some spikes, the overall pattern appears more random and less structured 
# than previous models.
# 
# The ACF plot shows dramatic improvement with virtually no significant 
# autocorrelation at any lag except lag 0. This indicates the temporal dependence 
# has been effectively captured by the model.
# 
# The PACF plot shows only scattered, small significant lags with no clear pattern. 



max_order <- 10
aic_values <- numeric(max_order+1)  
bic_values <- numeric(max_order+1)
rmse_values <- numeric(max_order+1)
model_list <- list()

for (p in 0:max_order) {
  if (p == 0) {
    # Order 0 model: just the current inputs, no autoregressive terms
    formula_str <- "Ph ~ Tdelta + Gv"
  } else {
    # Build formula string based on order p
    ph_terms <- paste0("Ph.l", 1:p, collapse = " + ")
    if (p > 0) {
      tdelta_lags <- paste0("Tdelta.l", 1:p-1, collapse = " + ")
      gv_lags <- paste0("Gv.l", 1:p-1, collapse = " + ")
      
      formula_str <- paste0("Ph ~ ", ph_terms," + ", tdelta_lags, " + ", gv_lags)
    } else {
      formula_str <- paste0("Ph ~ ", ph_terms, " + ", tdelta_current, " + ", gv_current)
    }
  }
  
  formula_obj <- as.formula(formula_str)
  cat("Fitting model of order", p, "with formula:", formula_str, "\n")
  
  # Fit the model
  model <- lm(formula_obj, data = train)
  model_list[[p+1]] <- model  # +1 because R lists are 1-indexed
  predictions <- predict(model, newdata = test)
  residuals <- test$Ph - predictions
  rmse_values[p+1] <- sqrt((1/64) * sum(residuals^2))
  # Store AIC and BIC values
  aic_values[p+1] <- AIC(model)
  bic_values[p+1] <- BIC(model)
  
  cat("Order", p, "- AIC:", aic_values[p+1], "BIC:", bic_values[p+1], "\n")
}

# Plot AIC and BIC values
par(mfrow = c(1, 1))
plot(x=0:max_order, y=aic_values, type="b", col="blue", pch=19, xlab="Model Order", ylab="AIC/BIC", main="AIC and BIC for Different Model Orders")
lines(x=0:max_order, y=bic_values, col="red", pch=19, type="b")
legend("topright", legend=c("AIC", "BIC"), col=c("blue", "red"), pch=19)


cat("Best model order based on AIC:", which.min(aic_values)-1, "\n")
cat("Best model order based on BIC:", which.min(bic_values)-1, "\n")

plot(x=0:max_order, y=rmse_values, type="b", col="black", pch=19, xlab="Model Order", ylab="RMSE", main="RMSE for Different Model Orders")

