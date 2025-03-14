library(astsa)
library(httpgd)
# Open server to view plots
hgd()

#1.4 Plot ACF up to k=30 for \phi_1 = -0.7 and \phi_2 = -0.2
x <- arima.sim(n=200, model=list(order=c(2,0,0), ar=c(-0.7, -0.2)))
acf(x, lag.max=30, lwd=2)

#2. Simulating seasonal processes 
n <- 300

plotit <- function(x){
  layout(rbind(1,2:3))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7,0))
  plot(x, ylab="X")
  acf(x, lag.max=50, lwd=2)
  pacf(x, lag.max=50, lwd=2)
}
#2.1 A (1,0,0) ×(0,0,0)12 model with the parameter φ1 = 0.6.
set.seed(11)
y <- arima.sim(n = n, model = list(ar = c(-0.6)))
plotit(y)


#2.2. A (0,0,0) ×(1,0,0)12 model with the parameter Φ1 = −0.9.
set.seed(12)
y <- arima.sim(n = 100, list(ar = c(rep(0,11), 0.9)))
plotit(y)

#2.3
#A (1,0,0) ×(0,0,1)12 model with the parameters φ1 = 0.9 and Θ1 = −0.7
set.seed(3)
y <- arima.sim(n = n, model = list(ar = c(-0.9), ma = c(rep(0,11), -0.7)))
plotit(y)

#2.4
#A (1,0,0) ×(1,0,0)12 model with the parameters φ1 = −0.6 and Φ1 = −0.8
set.seed(10)
y <- arima.sim(list(ar=c(0.6,rep(0,10),0.8,-0.6*0.8)),n)
plotit(y)

#2.5
#A (0,0,1) ×(0,0,1)12 model with the parameters θ1 = 0.4 and Θ1 = −0.8
set.seed(15)
y <- arima.sim(n = n, model = list(ma = c(0.4, rep(0,10), -0.8)))
plotit(y)


#2.6
#A (0,0,1) ×(1,0,0)12 model with the parameters θ1 = −0.4 and Φ1 = 0.7
set.seed(6)
y <- arima.sim(n = n, model = list(ma = c(-0.4), ar = c(rep(0,11), -0.7)))
plotit(y)

