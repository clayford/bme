# Biostats in epidemiology
# ch 10 work


# 10.1 --------------------------------------------------------------------

# Weibull hazard function
time <- seq(0,4,0.01)
weibull.hazard <- function(alpha, lambda, t){
  if(alpha <= 0 || lambda <= 0) stop("alpha and lambda must be positive.")
  alpha * lambda * (lambda * t)^(alpha - 1)
}

# wh1 <- weibull.hazard(alpha = 3.0, lambda = 0, t = time)
wh1 <- weibull.hazard(alpha = 3.0, lambda = 1, t = time)
wh2 <- weibull.hazard(alpha = 1.5, lambda = 1, t = time)
wh3 <- weibull.hazard(alpha = 1.0, lambda = 1, t = time)
wh4 <- weibull.hazard(alpha = 0.5, lambda = 1, t = time)

# Fig 10.1(a)
plot(time, wh1, ylim = c(0,4), ylab = "Hazard", xlab = "Time", type="l")
lines(time, wh2, lty = 2)
lines(time, wh3, lty = 3)
lines(time, wh4, lty = 4)
legend("topright", legend = paste("alpha =", c(3.0,1.5,1.0,0.5)), lty= 1:4)

# Weibull survival function
weibull.survival <- function(alpha, lambda, t){
  if(alpha <= 0 || lambda <= 0) stop("alpha and lambda must be positive.")
  exp(-(lambda * t)^alpha)
}

ws1 <- weibull.survival(alpha = 3.0, lambda = 1, t = time)
ws2 <- weibull.survival(alpha = 1.5, lambda = 1, t = time)
ws3 <- weibull.survival(alpha = 1.0, lambda = 1, t = time)
ws4 <- weibull.survival(alpha = 0.5, lambda = 1, t = time)
# Fig 10.1(b)
plot(time, ws1, ylim = c(0,1), ylab = "Survival Probability", xlab = "Time", type="l")
lines(time, ws2, lty = 2)
lines(time, ws3, lty = 3)
lines(time, ws4, lty = 4)
legend("topright", legend = paste("alpha =", c(3.0,1.5,1.0,0.5)), lty= 1:4)

# Example 10.1

#   survreg's scale  =    1/(rweibull shape)
#   survreg's intercept = log(rweibull scale)

dat <- read.csv("BreastCancer.csv")
library(flexsurv)
fitw <- flexsurvreg(formula = Surv(time, status) ~ 1, data = dat, dist="weibull")
# Fig 10.2(b)
plot(fitw, ci=FALSE)
fitws <- summary(fitw)
fitws[[1]]$est

fite <- flexsurvreg(formula = Surv(time, status) ~ 1, data = dat, dist="exp")
# Fig 10.2(a)
plot(fite, ci=FALSE)
fites <- summary(fite)
fites[[1]]$est

# this does not agree with book
library(exptest)
co.exp.test(x = fites[[1]]$est)


library(survival)
m1 <- survreg(Surv(time, status) ~ 1, data=dat, dist = "weibull")
sm1 <- summary(m1)
1/sm1$scale # alpha-hat (shape)
1/exp(sm1$coefficients[1]) # lambda-hat (rate)

# cox-oakes test of exponentiality
# Null: alpha = 1, Weibull shape parameter is 1, hazard function is constant
d <- sum(dat$status) # number of deaths in cohort
n <- sum(dat$time) # total amount of time cohort was under observation
lam0.hat <- d/n # deaths per person-month
si <- log(lam0.hat * dat$time)
ei.hat <- lam0.hat * dat$time

num <- (d + sum(si * (dat$status - ei.hat)))^2
denom <- d + sum((si^2)*ei.hat) - (sum(si*ei.hat))^2 / d
X.co <- num/denom
pchisq(X.co, df = 1, lower.tail = FALSE)
# provides moderate evidence that exponential assumption may not be satisfied

# Funciton: cox-oakes test
cox.oakes <- function(time,status) {
  dname <- deparse(substitute(time))
  d <- sum(status) # number of deaths in cohort
  n <- sum(time) # total amount of time cohort was under observation
  lam0.hat <- d/n # deaths per person-month
  si <- log(lam0.hat * time)
  ei.hat <- lam0.hat * time
  
  num <- (d + sum(si * (status - ei.hat)))^2
  denom <- d + sum((si^2)*ei.hat) - (sum(si*ei.hat))^2 / d
  X.co <- num/denom
  p.value <- pchisq(X.co, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = c(statistic = X.co), p.value = p.value, 
               method = "Cox-Oakes Test of Exponentiality", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

cox.oakes(time = dat$time, status = dat$status)
