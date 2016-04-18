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
  names(lam0.hat) <- "hazard rate"
  si <- log(lam0.hat * time)
  ei.hat <- lam0.hat * time
  
  num <- (d + sum(si * (status - ei.hat)))^2
  denom <- d + sum((si^2)*ei.hat) - (sum(si*ei.hat))^2 / d
  X.co <- num/denom
  p.value <- pchisq(X.co, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = c(statistic = X.co), p.value = p.value,
               estimate = lam0.hat,
               method = "Cox-Oakes Test of Exponentiality", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

cox.oakes(time = breast.survival$time, status = breast.survival$status)
c.out <- cox.oakes(time = breast.survival$time, status = breast.survival$status)
str(c.out)

# example 10.2
dat2 <- read.csv("ovarian.csv")

# Weibull model
m2 <- survreg(Surv(time, status) ~ 1, data=dat2, dist = "weibull", subset = grade=="High")
sm2 <- summary(m2)
1/sm2$scale # alpha-hat (shape)
1/exp(sm2$coefficients[1]) # lambda-hat (rate)

# Exponential model
m3 <- survreg(Surv(time, status) ~ 1, data=dat2, dist = "exponential", subset = grade=="High")
sm3 <- summary(m3)
1/exp(sm3$coefficients[1]) # lambda-hat (rate); alpha-hat (shape) is 1

# Fig 10.3(a)
library(flexsurv)
fite <- flexsurvreg(formula = Surv(time, status) ~ 1, data = dat2, dist="exp")
plot(fite, ci=FALSE)

# Fig 10.3(b)
fitw <- flexsurvreg(formula = Surv(time, status) ~ 1, data = dat2, dist="weibull")
plot(fitw, ci=FALSE)

cox.oakes(time = dat2$time[dat2$grade=="High"], status = dat2$status[dat2$grade=="High"])

# exact methods for single sample
# Example 10.3

# with n assumed constant (where n is total person time in study), d is a
# poisson random variable with parameter nu = lambda * n

# Null: lambda_0 = 0.4
d <- 2
n <- 10

ppois(q = 2, lambda = n*0.4) * 2

# Table 10.1
d <- 0:12
`D_=_d` <- dpois(d, lambda = n * 0.4)
`D_<=_d` <- ppois(d, lambda = n * 0.4)
`D_>=_d` <- ppois(d - 1, lambda = n * 0.4, lower.tail = FALSE)
cbind(d, `D_=_d`, `D_<=_d`, `D_>=_d`)


1 - ppois(q = 1, lambda = 10*0.024)
ppois(q = 2, lambda = 10*0.723)

f1 <- function(x)1 - ppois(q = 1, lambda = 10*x) - 0.025
f2 <- function(x)ppois(q = 2, lambda = 10*x) - 0.025
uniroot(f1, interval = c(0,10))
uniroot(f2, interval = c(0,10))

exact.rate.ci <- function(d,n,ci=0.95){
  alpha <- (1-ci)/2
  f1 <- function(x)1 - ppois(q = d-1, lambda = n*x) - alpha
  f2 <- function(x)ppois(q = d, lambda = n*x) - alpha
  f1.out <- uniroot(f1, interval = c(0,n))
  f2.out <- uniroot(f2, interval = c(0,n))
  list(lower = f1.out$root, upper = f2.out$root)
}
exact.rate.ci(d = 2, n = 10)
exact.rate.ci(d = 2, n = 10, ci=0.90)


exact.rate.test <- function(d, n, null=0, conf.level = 0.95){
  # estimate
  est <- d/n
  # CI
  alpha <- (1-conf.level)/2
  f1 <- function(x)1 - ppois(q = d-1, lambda = n*x) - alpha
  f2 <- function(x)ppois(q = d, lambda = n*x) - alpha
  f1.out <- uniroot(f1, interval = c(0,n))
  f2.out <- uniroot(f2, interval = c(0,n))
  # test
  p.value <- min(ppois(q = d, lambda = n*null) * 2, 1)
  list(estimate = est, null = null, p.value = round(p.value,3),
       conf.int = c(round(f1.out$root,3), round(f2.out$root,3)))
}

exact.rate.test(d = 2, n = 10, null = 0.4)

# asympotic methods for single sample
rate.test <- function(d, n, null=0, conf.level = 0.95, explicit=TRUE){
  # estimate
  est <- d/n
  # CI
  alpha <- (1-conf.level)/2
  if(explicit){
    ci <- est + c(-1,1)*((qnorm(1 - alpha)*sqrt(d))/n)  
  } else {
    a <- n^2
    b <- -n*(2*d + qnorm(1 - alpha)^2)
    c <- d^2
    ci <- (-b + c(-1,1)*sqrt(b^2 - 4*a*c))/(2*a)
  }
  # test
  test.statistic <- ((d - null*n)^2)/(null*n)
  p.value <- pchisq(test.statistic, df = 1, lower.tail = FALSE)
  list(estimate = est, null = null, test.statistic = round(test.statistic,3), 
       p.value = round(p.value,3),
       conf.int = round(ci,3))
}
# Examples 10.5 and 10.6
rate.test(d = 2, n = 10, null = 0.4)
rate.test(d = 5, n = 25, null = 0.4)
rate.test(d = 10, n = 50, null = 0.4)

rate.test(d = 2, n = 10, null = 0.4, explicit = FALSE)
rate.test(d = 5, n = 25, null = 0.4, explicit = FALSE)
rate.test(d = 10, n = 50, null = 0.4, explicit = FALSE)

#Example 10.7
exact.rate.test(d = sum(dat$status), n = sum(dat$time))
rate.test(d = sum(dat$status), n = sum(dat$time))
rate.test(d = sum(dat$status), n = sum(dat$time), explicit = FALSE)



# 10.2 --------------------------------------------------------------------

# Poisson methods for unstratified survival data

