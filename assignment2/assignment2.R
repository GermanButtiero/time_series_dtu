#2. Simulating seasonal processes 
library(astsa)

n <- 100

plotit <- function(x){
  layout(rbind(1,2:3))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7,0))
  plot(x, ylab="X")
  acf(x, lag.max=50, lwd=2)
  pacf(x, lag.max=50, lwd=2)
}
#2.1 A (1,0,0) ×(0,0,0)12 model with the parameter φ1 = 0.6.
set.seed(1)
y <- arima.sim(n = n, list(ar = c(0.6)))
plotit(y)

#2.2. A (0,0,0) ×(1,0,0)12 model with the parameter Φ1 = −0.9.
set.seed(2)
y <- sarima.sim(n = n, sar=c(-0.9), S=12)
plotit(y)

#2.3
#A (1,0,0) ×(0,0,1)12 model with the parameters φ1 = 0.9 and Θ1 = −0.7
set.seed(3)
y <- sarima.sim(n = n, ar= 0.9, sma = -0.7, S=12)
plotit(y)

#2.4
#A (1,0,0) ×(1,0,0)12 model with the parameters φ1 = −0.6 and Φ1 = −0.8
set.seed(4)
y <- sarima.sim(n = n, ar= -0.6, sar=c(-0.8), S=12)
plotit(y)

#2.5
#A (0,0,1) ×(0,0,1)12 model with the parameters θ1 = 0.4 and Θ1 = −0.8
set.seed(5)

y <- sarima.sim(n = n, ma= 0.4, sma = -0.8, S=12)
plotit(y)


#2.6
#A (0,0,1) ×(1,0,0)12 model with the parameters θ1 = −0.4 and Φ1 = 0.7
set.seed(6)
y <- sarima.sim(n = n, ma= -0.4, sar=c(0.7), S=12)
plotit(y)
